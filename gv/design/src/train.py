import os
os.environ['KMP_DUPLICATE_LIB_OK'] = 'TRUE'
import torch
import torch.optim as optim
from torch.utils.data import DataLoader
from data_loader import GvpDataset, collate_fn, VOCAB_SIZE
from models.vae_gvp import GVAE, vae_loss
import argparse

def train_vae(args):
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    os.makedirs('checkpoints', exist_ok=True)

    dataset = GvpDataset(args.fasta, args.msa, max_len=args.max_len, chain_type=args.chain_type)
    dataloader = DataLoader(dataset, batch_size=args.batch_size, 
                           shuffle=True, collate_fn=collate_fn)

    model = GVAE(vocab_size=VOCAB_SIZE,
                 num_species=dataset.num_species,
                 latent_dim=args.latent_dim,
                 target_len_A=args.target_len_A,
                 target_len_C=args.target_len_C).to(device)

    optimizer = optim.Adam(model.parameters(), lr=args.lr)

    for epoch in range(args.epochs):
        model.train()
        total_loss = 0
        num_batches = 0

        for batch in dataloader:
            x = batch['seq'].to(device)
            species = batch['species'].to(device)

            optimizer.zero_grad()
            outputs, mu, logvar = model(x, species)

            target = x[:, 1:]
            min_len = min(outputs.size(1), target.size(1))
            outputs = outputs[:, :min_len, :]
            target = target[:, :min_len]

            loss, recon, kl = vae_loss(outputs, target, mu, logvar, 
                                        kl_weight=args.kl_weight)

            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
            optimizer.step()

            total_loss += loss.item()
            num_batches += 1

        avg_loss = total_loss / max(num_batches, 1)
        print(f"Epoch {epoch}: Loss {avg_loss:.4f}")

        if epoch % 10 == 0:
            torch.save(model.state_dict(), f"checkpoints/vae_epoch{epoch}.pt")

    torch.save(model.state_dict(), "checkpoints/vae_final.pt")
    print("VAE训练完成！")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="VAE-only reference trainer; see gv/README.md before training")
    parser.add_argument('--model', choices=['vae'], default='vae')
    parser.add_argument('--fasta', required=True, help='Explicit GvpA training FASTA path')
    parser.add_argument('--msa', default=None, help='Optional MSA file or directory')
    parser.add_argument('--chain_type', nargs='+', default=None)
    parser.add_argument('--epochs', type=int, default=100)
    parser.add_argument('--batch_size', type=int, default=32)
    parser.add_argument('--lr', type=float, default=1e-4)
    parser.add_argument('--latent_dim', type=int, default=64)
    parser.add_argument('--kl_weight', type=float, default=0.1)
    parser.add_argument('--max_len', type=int, default=512)
    parser.add_argument('--target_len_A', type=int, default=100)
    parser.add_argument('--target_len_C', type=int, default=500)
    train_vae(parser.parse_args())
