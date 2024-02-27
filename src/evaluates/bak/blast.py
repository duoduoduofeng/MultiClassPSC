import ssl
from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

# 定义两个氨基酸序列
# sequence1 = "MNIRPLHDRVIVKRKEVETKSAGGIVLTGSAAAKSTRGEVLAVGNGRILENGEVKPLDVKVGDIVIFNDGYGVKSEKIDNEEVLIMSESDILAIVEA"
# sequence2 = "IRPLHDRVIVKRKEVETKSAGGIVLTGSAAAKSTRGEVLAVGNGRILENGEVKPLDVKVGDIVIFNDGYGVKSEKIDNEEVLIMSESDILAIV"

# sequence1 = "MNIRPLHDRVIVKRKEVETKSAGGIVLTGSAAAKSTRGEVLAVGNGRILENGEVKPLDVKVGDIVIFNDGYGVKSEKIDNEEVLIMSESDILAIVEA"
# sequence2 = "QQLPIRAVGEYVILVSEPELCVVHSVGPDVPEGFCEVGDLTSLPVGQIRNVPHPFVALGLKQPKEIKQKFVTCHYKAIPCLYK"

# class: 1
sequence1 = "MNIRPLHDRVIVKRKEVETKSAGGIVLTGSAAAKSTRGEVLAVGNGRILENGEVKPLDVKVGDIVIFNDGYGVKSEKIDNEEVLIMSESDILAIVEA"
sequence2 = "MIKPLGDRVVVKRIVLPDTAKEKPQKGKVIAVGTGRVLENGQRVPLEVKEGDIVVFAKYGGTEIEIDGEEYVILSERDLLAVLQ"

# 禁用SSL证书验证
# ssl._create_default_https_context = ssl._create_unverified_context

# # 执行BLAST搜索
# result_handle = NCBIWWW.qblast("blastp", "nr", sequence1)

# # 解析BLAST结果
# blast_record = NCBIXML.read(result_handle)

# # 获取第一个比对结果的比对分数
# alignment = blast_record.alignments[0]
# print(alignment)
# hsp = alignment.hsps[0]
# print("Alignment Score:", hsp.score)

# 构建序列字典
sequences = {
    "sequence1": sequence1,
    "sequence2": sequence2
}

# 执行 Blastp 搜索，同时传递两个序列
# result_handle = NCBIWWW.qblast("blastp", "nr", sequences)

# 解析 Blastp 结果
# blast_record = NCBIXML.read(result_handle)

# # Convert to a dictionary
# blast_dict = blast_record.to_dict()


# 获取第一个比对结果的比对信息
# alignment = blast_record.alignments[0]
# hsp = alignment.hsps[0]

# # 计算每个比对的身份值
# alignment_length = hsp.align_length
# identities = hsp.identities
# per_identity = (identities / alignment_length) * 100

# print("Per Identity:", per_identity)

# 获取每个比对的 per identity 值
# alignment_count = 0
# hsp_count = 0

# for alignment in blast_record.alignments:
# 	alignment_count += 1
# 	hsp_count = 0
# 	for hsp in alignment.hsps:
# 		hsp_count += 1
# 		per_identity = hsp.identities / hsp.align_length
# 		print(f"alignment_count = {alignment_count}, hsp_count = {hsp_count}, Per Identity: {per_identity}")

# 		alignment_score = hsp.score
# 		print(f"alignment_count = {alignment_count}, hsp_count = {hsp_count}, Alignment Score: {alignment_score}")

# # 存储所有比对的身份值和比对长度
# identities = []
# align_lengths = []

# # 获取所有比对的身份值和比对长度
# for alignment in blast_record.alignments:
# 	for hsp in alignment.hsps:
# 		alignment_length = hsp.align_length
# 		identity = (hsp.identities / alignment_length) * 100
# 		identities.append(identity)
# 		align_lengths.append(alignment_length)

# # 计算加权平均的 per identity
# if identities:
# 	total_identity = sum(identity * align_length for identity, align_length in zip(identities, align_lengths))
# 	total_align_length = sum(align_lengths)
# 	average_per_identity = total_identity / total_align_length
# 	print("Average Per Identity:", average_per_identity)
# else:
# 	print("No alignment results found")

# 计算相似度
# 这里简单地使用比对分数作为相似度指标，您可以根据需要选择其他指标
# 通常情况下，相似度可以使用比对分数、比对长度、或标准化后的比对分数来表示
# alignment_score = hsp.score
# similarity = alignment_score / max(len(sequence1), len(sequence2))

# print("Alignment Score:", alignment_score)
# print("Similarity:", similarity)


# 定义两个蛋白质序列
query_sequence = sequence1
subject_sequence = sequence2

# 将 subject_sequence 保存到一个 FASTA 文件中
with open("subject_sequence.fasta", "w") as f:
	f.write(">subject\n" + subject_sequence)

# 执行 BLAST 搜索，将一个序列作为查询序列，另一个作为数据库
result_handle = NCBIWWW.qblast("blastp", "nr", query_sequence, 
	db="subject_sequence.fasta", alignments=1, descriptions=1)

# 解析 BLAST 结果
blast_record = NCBIXML.read(result_handle)

# 获取第一个比对结果的比对信息
alignment = blast_record.alignments[0]
hsp = alignment.hsps[0]
alignment_length = hsp.align_length
identities = hsp.identities
# 计算相似度
similarity = (identities / alignment_length) * 100

print("Similarity:", similarity)


# 关闭结果句柄
result_handle.close()
