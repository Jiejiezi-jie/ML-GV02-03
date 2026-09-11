import torch
import os
from model import TransformerLM
from utils import BOS_ID, EOS_ID, decode, top_k_sampling, VOCAB

# 1. 硬件配置
device = "cuda" if torch.cuda.is_available() else "cpu"

# 2. 路径配置（匹配你的 AutoDL 目录截图）
# 确保与 train.py 保存的路径一致
CKPT_PATH = "gvp_ckpt/model.pt"  
# 最终产物 fasta 文件的保存路径
OUTPUT_PATH = "gvp_output/generated_gvp.fasta"  

# 确保输出文件夹存在，防止报错
os.makedirs("gvp_output", exist_ok=True)

# 3. 加载模型
print(f"正在从 {CKPT_PATH} 加载模型权重...")
model = TransformerLM(len(VOCAB))

# 使用 map_location 增加兼容性
try:
# 确保 generate.py 也能找到这个新位置
    model.load_state_dict(torch.load("/root/autodl-tmp/gvp_ckpt/model.pt", map_location=device))
except FileNotFoundError:
    print(f"错误：在 {CKPT_PATH} 未找到权重文件，请确认训练是否完成并保存。")
    exit()

model.to(device)
model.eval() # 切换到评估模式

def generate(max_len=512, temperature=0.8):
    """
    max_len: 建议设为 512 或更高，以兼容 GvpN 等长蛋白
    temperature: 采样温度，越高随机性越强
    """
    # 从起始符 <BOS> 开始生成
    seq = torch.tensor([[BOS_ID]], device=device)
    
    with torch.no_grad(): # 推理模式不计算梯度，节省显存
        for _ in range(max_len):
            # 获取模型对当前序列最后一个位置的预测 logits
            logits = model(seq)[:, -1, :] / temperature
            
            # 使用 top_k 采样选择下一个氨基酸 ID
            next_id = top_k_sampling(logits.squeeze(), k=10)
            
            # 将新生成的 ID 拼接到序列中
            seq = torch.cat([seq, torch.tensor([[next_id]], device=device)], dim=1)
            
            # 如果生成了结束符 <EOS>，则停止生成
            if next_id == EOS_ID:
                break
                
    # 将 ID 序列解码回氨基酸字母序列
    return decode(seq.squeeze().tolist())

# 4. 执行批量生成
print(f"开始批量生成 200 条 Gvp 序列...")
with open(OUTPUT_PATH, "w") as f:
    for i in range(200):
        seq = generate()
        # 写入标准 FASTA 格式
        f.write(f">gen_{i}\n{seq}\n")
        
        if (i + 1) % 50 == 0:
            print(f"已完成: {i+1}/200")

print(f"✅ 生成完毕！产物已存至: {OUTPUT_PATH}")