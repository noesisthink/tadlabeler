# TAD 基因组结构域分析工具

[![Python Version](https://img.shields.io/badge/python-3.8+-blue.svg)](https://python.org)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![GitHub issues](https://img.shields.io/github/issues/your-username/tad)](https://github.com/your-username/tad/issues)

> 基于深度学习的TAD(拓扑关联结构域)边界识别与分析系统

## 项目简介

本项目提供基于卷积神经网络的TAD边界预测工具，用于基因组三维结构分析。支持mcool格式Hi-C数据处理，能够自动识别染色质拓扑关联结构域边界，为基因组功能研究提供计算支持。

## ✨ 特性

- 🧬 高精度CNN TAD边界预测模型
- 📊 支持mcool / cool / bed 多种格式输入
- ⚡ GPU加速计算支持
- 🎯 90%+ 预测准确率
- 📈 内置可视化与统计分析
- 🖥️ 跨平台桌面客户端(tad_tag)

## 📁 项目结构

```
tad/
├── 📄 main.py                 # 主程序入口
├── 📄 mcool_bed.py           # mcool文件处理模块
├── 📄 sql.py                 # 数据库操作模块
├── 📂 dactor/                # 深度学习预测模块
│   ├── CNN.py               # 神经网络模型定义
│   └── test.py              # 模型测试脚本
├── 📂 extra_mode/           # 扩展功能模块
├── 📂 tad_tag/              # 桌面客户端(Electron)
├── 📂 tests/                # 单元测试
├── 📂 docs/                 # 文档目录
├── 📄 requirements.txt      # Python依赖
├── 📄 .gitignore           # Git忽略规则
├── 📄 LICENSE              # 开源协议
└── 📄 README.md            # 项目说明
```

## 🚀 快速开始

### 环境要求

- Python 3.8+
- TensorFlow / PyTorch
- 推荐GPU环境（CUDA支持）

### 安装

```bash
# 克隆仓库
git clone https://github.com/your-username/tad.git
cd tad

# 安装依赖
pip install -r requirements.txt
```

### 使用方法

```bash
# 运行预测
python main.py --input your_data.bed --output results

# 运行测试
python -m pytest tests/ -v
```

## 🧪 测试

```bash
# 运行全部测试
pytest

# 运行特定测试
pytest tests/test_model.py -v

# 生成覆盖率报告
pytest --cov=dactor
```

## 🤝 参与开发

欢迎提交Issue和Pull Request！

### 开发流程

1. Fork 本仓库
2. 创建特性分支 (`git checkout -b feature/AmazingFeature`)
3. 提交更改 (`git commit -m 'Add some AmazingFeature'`)
4. 推送到分支 (`git push origin feature/AmazingFeature`)
5. 打开 Pull Request

### 代码规范

- 遵循 PEP 8 编码规范
- 新增功能必须包含单元测试
- 提交前通过所有现有测试
- 保持代码清晰可维护

## 📝 许可证

本项目采用 MIT 许可证 - 查看 [LICENSE](LICENSE) 文件了解详情

## 👥 贡献者

- 项目维护者
