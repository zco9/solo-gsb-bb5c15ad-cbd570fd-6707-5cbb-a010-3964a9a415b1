# 化学方程式精确配平 API

基于 FastAPI 的化学方程式配平服务。服务端为反应物/生成物构建**元素守恒矩阵（含电荷行）**，用**有理数高斯消元**求零空间；仅当零空间为一维且基向量可化为全正时，返回最大公约数约简后的**最小正整数系数**及逐元素、电荷核算，否则返回 `422` 并说明原因（零空间维数或无法全正）。

另提供**化学计量接口**：在配平结果基础上，按内置原子量表（H、C、N、O、Na、Mg、Al、S、Cl、K、Ca、Fe、Cu、Mn）换算投料物质的量，以反应进度最小值判定（可并列的）限制反应物，返回理论反应进度、各生成物理论产量与各反应物剩余量。全程 `Decimal` 计算，输出按请求 `precision` 四舍五入。

还提供**经验式反推接口**：由元素质量百分比按原子量换算摩尔比，在 1～12 倍尺度内枚举互质正整数下标作为经验式候选，反算理论质量百分比择优；可选结合实测摩尔质量确定分子式。

## 快速启动（单条命令）

```bash
docker compose up --build
```

依赖随容器启动自动安装（`pip install -r requirements.txt` 后再拉起 uvicorn）。

- 服务地址：http://localhost:8000
- 交互式文档（Swagger UI）：http://localhost:8000/docs
- 健康检查：http://localhost:8000/health

## 配置

| 环境变量 | 默认值 | 说明 |
| --- | --- | --- |
| `PORT` | `8000` | 容器内服务监听端口 |

宿主机端口映射在 `docker-compose.yml` 的 `ports` 中调整（默认 `8000:8000`）。

## 分子式语法

| 语法 | 示例 |
| --- | --- |
| 元素与下标 | `H2O`、`Fe2O3` |
| 嵌套括号（`()` `[]` `{}` 可混用） | `Ca(OH)2`、`K4[Fe(CN)6]` |
| 点号结晶水（`.` `·` `*`） | `CuSO4.5H2O`、`Na2CO3·10H2O` |
| 插入符电荷（推荐，无歧义） | `SO4^2-`、`Fe^3+`、`NH4^+` |
| 纯符号尾电荷（每个 `+`/`-` 计 1） | `Na+`、`Cl-`、`Ca++`、`SO4--` |
| 电子（半反应） | `e-` |

> 注意：不带 `^` 时，符号前的数字属于前一个元素/基团——`NH4+` 解析为 NH₄⁺（电荷 +1）；要表示 Ca²⁺ 请写 `Ca^2+`。

## API

### `POST /api/v1/balance` — 配平

```bash
curl -s http://localhost:8000/api/v1/balance \
  -H 'Content-Type: application/json' \
  -d '{"reactants": ["KMnO4", "HCl"], "products": ["KCl", "MnCl2", "H2O", "Cl2"]}'
```

响应（`200`）：

```json
{
  "equation": "2KMnO4 + 16HCl -> 2KCl + 2MnCl2 + 8H2O + 5Cl2",
  "coefficients": [2, 16, 2, 2, 8, 5],
  "reactant_coefficients": [2, 16],
  "product_coefficients": [2, 2, 8, 5],
  "elements": [
    {"element": "Cl", "reactants": 16, "products": 16, "balanced": true},
    {"element": "H",  "reactants": 16, "products": 16, "balanced": true},
    {"element": "K",  "reactants": 2,  "products": 2,  "balanced": true},
    {"element": "Mn", "reactants": 2,  "products": 2,  "balanced": true},
    {"element": "O",  "reactants": 8,  "products": 8,  "balanced": true}
  ],
  "charge": {"reactants": 0, "products": 0, "balanced": true},
  "balanced": true
}
```

离子方程式（含电荷行守恒）：

```bash
curl -s http://localhost:8000/api/v1/balance \
  -H 'Content-Type: application/json' \
  -d '{"reactants": ["MnO4-", "H+", "e-"], "products": ["Mn^2+", "H2O"]}'
# -> MnO4- + 8H+ + 5e- -> Mn^2+ + 4H2O
```

无法配平时返回 `422`，并说明零空间维数或无法全正：

```json
{
  "detail": {
    "message": "零空间维数为 2（必须恰好为 1 才能唯一配平）",
    "null_space_dim": 2
  }
}
```

化学式语法错误返回 `400`。

### `POST /api/v1/verify` — 校验系数

```bash
curl -s http://localhost:8000/api/v1/verify \
  -H 'Content-Type: application/json' \
  -d '{
        "reactants": ["H2", "O2"],
        "products": ["H2O"],
        "reactant_coefficients": [2, 1],
        "product_coefficients": [2]
      }'
```

响应包含逐元素与电荷核算、`balanced` 结论，以及系数的 `gcd` 与是否已约简（`is_reduced`）。

### `POST /api/v1/stoichiometry` — 化学计量

在配平结果基础上，根据各反应物投料（单位仅支持 `mol` 或 `g`）计算限制反应物、理论反应进度、各生成物理论产量与各反应物剩余量。

- 内置原子量表（g/mol）：H 1.008、C 12.011、N 14.007、O 15.999、Na 22.98977、Mg 24.305、Al 26.982、S 32.06、Cl 35.45、K 39.0983、Ca 40.078、Fe 55.845、Cu 63.546、Mn 54.938（离子电荷的电子质量忽略不计）。
- 每种反应物必须恰好投料一次；以反应进度 ξ = n/ν 的最小值判定限制反应物（可并列多个）。
- 全程 `Decimal` 计算，所有小数字段按 `precision`（2～8，默认 4）四舍五入为**定长小数字符串**返回。

```bash
curl -s http://localhost:8000/api/v1/stoichiometry \
  -H 'Content-Type: application/json' \
  -d '{
        "reactants": ["H2", "O2"],
        "products": ["H2O"],
        "feed": [
          {"formula": "H2", "amount": 10,  "unit": "mol"},
          {"formula": "O2", "amount": 100, "unit": "g"}
        ],
        "precision": 4
      }'
```

响应（`200`）：

```json
{
  "equation": "2H2 + O2 -> 2H2O",
  "coefficients": [2, 1, 2],
  "reactant_coefficients": [2, 1],
  "product_coefficients": [2],
  "precision": 4,
  "extent": "3.1252",
  "limiting_reactants": ["O2"],
  "products": [
    {"formula": "H2O", "coefficient": 2, "theoretical_moles": "6.2504", "theoretical_grams": "112.6008"}
  ],
  "reactants": [
    {"formula": "H2", "coefficient": 2, "feed_moles": "10.0000", "feed_grams": "20.1600",
     "consumed_moles": "6.2504", "remaining_moles": "3.7496", "remaining_grams": "7.5592", "limiting": false},
    {"formula": "O2", "coefficient": 1, "feed_moles": "3.1252", "feed_grams": "100.0000",
     "consumed_moles": "3.1252", "remaining_moles": "0.0000", "remaining_grams": "0.0000", "limiting": true}
  ]
}
```

错误返回：

| 状态码 | 场景 |
| --- | --- |
| `400` | 化学式无法解析；缺失投料；同一物质重复投料；投料量非正数；投料物质不是反应物；元素无内置原子量 |
| `422` | 方程无法配平（零空间维数不为 1 或无法全正）；单位不是 `mol`/`g`；`precision` 超出 2～8 |

### `POST /api/v1/formula/infer` — 反推经验式 / 分子式

由元素**质量百分比**反推经验式，提供实测摩尔质量时进一步给出分子式。

- 算法：各元素百分比除以内置原子量得摩尔比，以最小摩尔值归一；在 **1～12 倍尺度**内四舍五入为正整数下标，仅保留**互质**（gcd = 1）下标作为经验式候选。
- 候选反算各元素理论质量百分比，要求**最大绝对偏差**（百分点）不超过请求 `tolerance`（默认 `0.1`，必须为正）；偏差相同时依次取**总原子数较少**、**Hill 顺序字符串较小**者；无合格候选返回 `422`。
- 输入百分比均须为**正数**，且总和与 `100` 的绝对偏差不超过 `tolerance`，否则返回 `400`；元素须在内置原子量表内。
- 提供 `molar_mass` 时，取 `实测值 / 经验式摩尔质量` **最接近的正整数倍** n 确定分子式（理论值 = n × 经验式摩尔质量），其与实测值的**相对偏差不得超过 `tolerance`%**，否则 `422`。
- 化学式按 **Hill 系统**输出（含碳时 C、H 在前，其余字母序；不含碳时全部字母序）；所有小数字段四舍五入为 6 位小数字符串。

仅反推经验式（如 Fe 69.94% / O 30.06%，摩尔比约 1:1.5，需放大到尺度 2）：

```bash
curl -s http://localhost:8000/api/v1/formula/infer \
  -H 'Content-Type: application/json' \
  -d '{"elements": {"Fe": "69.94", "O": "30.06"}, "tolerance": 0.1}'
```

响应（`200`）：

```json
{
  "empirical_formula": "Fe2O3",
  "subscripts": {"Fe": 2, "O": 3},
  "formula_scale": 2,
  "total_atoms": 5,
  "elements": [
    {"element": "Fe", "subscript": 2, "given_percent": "69.940000",
     "theoretical_percent": "69.943076", "deviation": "0.003076"},
    {"element": "O", "subscript": 3, "given_percent": "30.060000",
     "theoretical_percent": "30.056924", "deviation": "-0.003076"}
  ],
  "max_abs_deviation": "0.003076",
  "empirical_molar_mass": "159.687000",
  "measured_molar_mass": null,
  "molecular_multiplier": null,
  "molecular_formula": null,
  "molecular_molar_mass": null,
  "molar_mass_relative_deviation_percent": null
}
```

结合实测摩尔质量确定分子式（葡萄糖：经验式 CH₂O，M = 30.026，实测 180.16 → n = 6）：

```bash
curl -s http://localhost:8000/api/v1/formula/infer \
  -H 'Content-Type: application/json' \
  -d '{
        "elements": {"H": "6.71", "C": "40.00", "O": "53.29"},
        "tolerance": 0.1,
        "molar_mass": 180.16
      }'
```

响应（`200`，节选）：

```json
{
  "empirical_formula": "CH2O",
  "subscripts": {"C": 1, "H": 2, "O": 1},
  "formula_scale": 1,
  "max_abs_deviation": "0.006179",
  "empirical_molar_mass": "30.026000",
  "measured_molar_mass": "180.160000",
  "molecular_multiplier": 6,
  "molecular_formula": "C6H12O6",
  "molecular_molar_mass": "180.156000",
  "molar_mass_relative_deviation_percent": "0.002220"
}
```

错误返回：

| 状态码 | 场景 |
| --- | --- |
| `400` | `elements` 为空；存在表外元素；百分比非正或非有限数；百分比总和与 100 的偏差超过 `tolerance`；`tolerance`/`molar_mass` 非正 |
| `422` | 1～12 尺度内无最大绝对偏差不超过 `tolerance` 的候选；提供摩尔质量时最近正整数倍的相对偏差超过 `tolerance`%；字段类型/缺失等模型校验错误 |

## 本地开发

```bash
pip install -r requirements-dev.txt
uvicorn app.main:app --reload --port 8000
pytest                 # 运行测试
```

## 项目结构

```
app/
  parser.py         # 分子式解析：嵌套括号、结晶水、离子电荷
  balancer.py       # 守恒矩阵（含电荷行）、有理数 RREF、零空间、整数化
  stoichiometry.py  # 内置原子量表、投料换算、限制反应物、产量与剩余量（Decimal）
  formula_infer.py  # 质量百分比反推经验式、摩尔质量定分子式（Decimal）
  schemas.py        # 请求/响应模型
  main.py           # FastAPI 路由
tests/test_api.py / tests/test_stoichiometry.py / tests/test_formula_infer.py
Dockerfile / docker-compose.yml
```
