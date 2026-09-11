import torch
from torch.utils.data import DataLoader
from dataset import FastaDataset, collate_fn
from model import TransformerLM
from utils import PAD_ID, VOCAB
from tqdm import tqdm

device = "cuda" if torch.cuda.is_available() else "cpu"

dataset = FastaDataset("/root/autodl-tmp/gvp_data/gvp_clean.fasta")
loader = DataLoader(dataset, batch_size=32, shuffle=True, collate_fn=collate_fn)

model = TransformerLM(len(VOCAB)).to(device)
optimizer = torch.optim.AdamW(model.parameters(), lr=3e-4)
criterion = torch.nn.CrossEntropyLoss(ignore_index=PAD_ID)

for epoch in range(40):
    model.train()
    total_loss = 0

    for batch in tqdm(loader, desc=f"Epoch {epoch}"):
        batch = batch.to(device)
        logits = model(batch[:, :-1])
        loss = criterion(
            logits.reshape(-1, logits.size(-1)),
            batch[:, 1:].reshape(-1)
        )

        optimizer.zero_grad()
        loss.backward()
        optimizer.step()

        total_loss += loss.item()

    print(f"Epoch {epoch} | Loss: {total_loss / len(loader):.4f}")

torch.save(model.state_dict(), "/root/autodl-tmp/gvp_ckpt/model.pt")
