# FBcount
Tools for Feature Barcode Counting. The Lite version of FBranger contributed by collaborator Fengrui Sui.


## Install

To access the latest version of the code, clone and install the repository:

```bash
git clone git@github.com:yuanfu-bio/FBcount.git
```
For easiest use, [create a conda environment](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-with-commands) and install requirements:

```bash
conda create -n <env_name> python=3.11
conda install -c conda-forge jq
conda activate <env_name>
cd FBcount
pip install -r ./requirements.txt
```

## Usage


## 📌 To-Do List
### 🐛 Features
1. 优化饱和度计算
- [x]  所有barcode统一抽样，不再区分barcode，避免极微量情况下的UMI在不同barcode间的抽样异常
- [x]  添加饱和度抽样报告（出图）

### 🛠️ Bug Fixes
- [x]  生成抽样校正后的最终count矩阵

### 📖 Documentation
- [x] Write user guide for installation
- [ ] Write usage

### 🚀 Optimization
1.  优化代码结构
- [x]  整理通用函数
- [x]  完善日志信息
- [ ]  优化json输入（包括白名单及FB information）
- [ ]  统一代码风格
2. 代码环境配置优化
- [x]  生成requirements.txt
- [x]  整理本地环境shebang