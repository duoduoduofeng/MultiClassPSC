from Bio import pairwise2
from Bio.Align import substitution_matrices


# 定义两个氨基酸序列
seq1 = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGAAQQQVVGSLAQALDWGK"
seq2 = "MKVLYLAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGAAQQQVVGSLAQALDWKK"

# 获取BLOSUM62打分矩阵
# matrix = matlist.blosum62
matrix = substitution_matrices.load('BLOSUM62')

# 全局对齐 - Needleman-Wunsch算法
global_alignment = pairwise2.align.globalds(seq1, seq2, matrix, -10, -0.5)[0]  # gap开启和扩展惩罚

# 局部对齐 - Smith-Waterman算法
local_alignment = pairwise2.align.localds(seq1, seq2, matrix, -10, -0.5)[0]

# 计算相似度值
def calculate_similarity(alignment, length):
    score, max_score = alignment.score, length * matrix[('A', 'A')]  # 假设最大得分为序列长度乘以最高单个匹配得分
    similarity = score / max_score
    return similarity

# 计算全局和局部对齐的相似度
global_similarity = calculate_similarity(global_alignment, max(len(seq1), len(seq2)))
local_similarity = calculate_similarity(local_alignment, max(len(seq1), len(seq2)))

print(global_similarity)
print(local_similarity)
