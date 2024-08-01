# MOLcreat

CORE | AAAI2020：分子自动优化模型
https://blog.csdn.net/u012325865/article/details/104489308

目标导向分子图生成的图卷积策略网络
https://zhuanlan.zhihu.com/p/99412778

PN-18: GCPN for Molecular Generation (NIPS 2018)
https://zhuanlan.zhihu.com/p/70134702?from_voters_page=true

DGL | 基于深度图学习框架DGL的分子图生成库
https://docs.dgl.ai/tutorials/models/1_gnn/9_gat.html

斯坦福教授ICLR演讲：图网络最新进展GraphRNN和GCPN（附PPT下载）
https://baijiahao.baidu.com/s?id=1635036452161379749&wfr=spider&for=pc

Nat. Mach. Intell. | 利用条件循环神经网络生成特定性质分子
https://blog.csdn.net/u012325865/article/details/106309160/?utm_medium=distribute.pc_relevant.none-task-blog-title-1&spm=1001.2101.3001.4242

MoFlow:AnInvertibleFlowModelforGenerating MolecularGraphs 
https://arxiv.org/pdf/2006.10137.pdf

https://arxiv.org/pdf/1811.12823.pdf
MolecularSets(MOSES):ABenchmarkingPlatform forMolecularGenerationModels

JunctionTreeVariationalAutoencoderforMolecularGraphGeneration
http://proceedings.mlr.press/v80/jin18a/jin18a.pdf

GraphNVP:AnInvertibleFlowModelfor GeneratingMolecularGraphs
https://arxiv.org/pdf/1905.11600.pdf

https://openreview.net/group?id=ICLR.cc/2021/Conference


```

from rdkit import Chem
from rdkit.Chem import AllChem
 
# 定义一个生成分子的函数
def generate_molecule(smiles_list, max_attempts=1000):
    mol = Chem.RWMol()
    for smiles in smiles_list:
        # 将SMILES片段转换为分子
        frag = Chem.MolFromSmiles(smiles)
        if frag is None:
            raise ValueError(f"无法从SMILES片段{smiles}生成分子")
        # 尝试添加分子片段到构造的分子上
        added = False
        for _ in range(max_attempts):
            if mol.addMol(frag):
                added = True
                break
        if not added:
            raise ValueError(f"无法添加SMILES片段{smiles}到分子构造中")
    # 尝试生成一个可能的分子结构
    AllChem.EmbedMultipleConfs(mol.GetMol(), numConfs=1, randomSeed=0xDEADBEEF)
    return mol
 
# 使用示例
smiles_list = ['C', 'CC', 'O']
molecule = generate_molecule(smiles_list)
# 打印分子的SMILES表示
print(Chem.MolToSmiles(molecule.GetMol()))

这段代码首先定义了一个generate_molecule函数，它接受一个SMILES字符串列表，然后尝试将每个SMILES片段添加到一个分子对象上。如果添加成功，它会尝试为这个分子生成一个结构。这个过程使用了一个循环和错误处理，以确保即使在面对无法添加的情况下，也不会导致程序崩溃。最后，它打印出生成的分子的SMILES表示。

```


