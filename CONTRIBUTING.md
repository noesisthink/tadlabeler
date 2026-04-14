# 贡献指南

欢迎各位同学参与本项目开发！请先阅读以下规范，让我们的协作更顺畅。

## 📋 开发规范

### 代码风格
- 遵循 **PEP 8** Python编码规范
- 4空格缩进，禁止使用Tab
- 函数和类必须有文档注释
- 复杂逻辑必须添加注释说明

### Git提交规范
```
feat: 新增功能
fix: 修复bug
docs: 文档更新
style: 代码格式修改
refactor: 代码重构
test: 测试相关
chore: 构建/工具变动
```

### 分支管理
- `main` 主分支（稳定代码）
- `dev` 开发分支
- 每个人在自己的 `feature/xxx` 分支开发
- 开发完成后提交Pull Request到dev分支

## 🚩 开发流程

1. 从dev分支创建自己的功能分支
```bash
git checkout dev
git checkout -b feature/your-feature-name
```

2. 开发完成后运行测试
```bash
pytest
```

3. 提交代码
```bash
git add .
git commit -m "feat: 你的提交说明"
git push origin feature/your-feature-name
```

4. 在Github上创建Pull Request

## ✅ Pull Request要求
- 必须通过所有测试
- 新增功能必须附带单元测试
- 代码不能有语法错误
- 必须描述清楚改动内容

## 📞 遇到问题
- 先查看Issue列表是否已有相关问题
- 没有的话新建Issue详细描述问题
- 或者在群里讨论

感谢大家的贡献！🎉