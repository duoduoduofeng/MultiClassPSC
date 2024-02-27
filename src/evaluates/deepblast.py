import requests

def deepblast(sequence1, sequence2):
    # 设置 DeepBLAST API 的 URL
    url = "https://api.deepblast.org/api/v1/align"
    
    # 准备 POST 请求的数据
    data = {
        "sequence_a": sequence1,
        "sequence_b": sequence2
    }

    # 发送 POST 请求
    response = requests.post(url, json=data)
    
    # 解析响应并获取相似度分数
    if response.status_code == 200:
        result = response.json()
        similarity_score = result["score"]
        return similarity_score
    else:
        print("请求失败:", response.status_code)
        return None

# 示例氨基酸序列
sequence1 = "MVHQAQLSYNYDKYPLEPHDLNLLRDTKVPFSFEVSPEVTVMKPNVQLQGSDPYLVGIDSVPTKVDVLEMSGVTTTKVDEIVPFKGSQFTADPSEKLKGISFPLSQETFSQVVTSS"
sequence2 = "MVHQAVLSYNQDKYPLEPHDFQVSGEVTVMKPNVQLQGSDPYLVGIDSVPTKVDVLEMSGVTTTKVDEIVPFKGSQFTADPSEKLKGISFPLSQETFSQVVTSS"

# 使用 DeepBLAST 比对两个序列并计算相似度
similarity_score = deepblast(sequence1, sequence2)
if similarity_score is not None:
    print("相似度分数:", similarity_score)
