"""请求 / 响应模型。"""
from __future__ import annotations

from decimal import Decimal
from typing import Literal

from pydantic import BaseModel, Field


class BalanceRequest(BaseModel):
    reactants: list[str] = Field(
        min_length=1,
        description="反应物分子式列表",
        examples=[["KMnO4", "HCl"]],
    )
    products: list[str] = Field(
        min_length=1,
        description="生成物分子式列表",
        examples=[["KCl", "MnCl2", "H2O", "Cl2"]],
    )


class ElementAccount(BaseModel):
    element: str
    reactants: int = Field(description="反应物侧该元素原子总数")
    products: int = Field(description="生成物侧该元素原子总数")
    balanced: bool


class ChargeAccount(BaseModel):
    reactants: int = Field(description="反应物侧总电荷")
    products: int = Field(description="生成物侧总电荷")
    balanced: bool


class BalanceResponse(BaseModel):
    equation: str
    coefficients: list[int] = Field(description="全部系数（先反应物后生成物）")
    reactant_coefficients: list[int]
    product_coefficients: list[int]
    elements: list[ElementAccount] = Field(description="逐元素核算")
    charge: ChargeAccount = Field(description="电荷核算")
    balanced: bool


class VerifyRequest(BaseModel):
    reactants: list[str] = Field(min_length=1, examples=[["H2", "O2"]])
    products: list[str] = Field(min_length=1, examples=[["H2O"]])
    reactant_coefficients: list[int] = Field(min_length=1, examples=[[2, 1]])
    product_coefficients: list[int] = Field(min_length=1, examples=[[2]])


class VerifyResponse(BaseModel):
    balanced: bool
    equation: str
    elements: list[ElementAccount]
    charge: ChargeAccount
    gcd: int = Field(description="系数的最大公约数")
    is_reduced: bool = Field(description="系数是否已约简（gcd == 1）")
    message: str


class FeedItem(BaseModel):
    formula: str = Field(
        description="投料物质分子式（必须是反应物之一）", examples=["H2"]
    )
    amount: Decimal = Field(description="投料量，必须为正数", examples=["2.5"])
    unit: Literal["mol", "g"] = Field(description="投料单位，仅支持 mol 或 g")


class StoichiometryRequest(BaseModel):
    reactants: list[str] = Field(
        min_length=1,
        description="反应物分子式列表",
        examples=[["H2", "O2"]],
    )
    products: list[str] = Field(
        min_length=1,
        description="生成物分子式列表",
        examples=[["H2O"]],
    )
    feed: list[FeedItem] = Field(
        description="各反应物的投料（每种反应物恰好一条，单位 mol 或 g）"
    )
    precision: int = Field(
        default=4,
        ge=2,
        le=8,
        description="输出小数位数（四舍五入），2～8，默认 4",
    )


class ProductYield(BaseModel):
    formula: str
    coefficient: int = Field(description="配平系数")
    theoretical_moles: str = Field(description="理论产量（mol，定精度小数字符串）")
    theoretical_grams: str = Field(description="理论产量（g，定精度小数字符串）")


class ReactantRemainder(BaseModel):
    formula: str
    coefficient: int = Field(description="配平系数")
    feed_moles: str = Field(description="投料量（mol）")
    feed_grams: str = Field(description="投料量（g）")
    consumed_moles: str = Field(description="消耗量（mol）")
    remaining_moles: str = Field(description="剩余量（mol）")
    remaining_grams: str = Field(description="剩余量（g）")
    limiting: bool = Field(description="是否为限制反应物")


class StoichiometryResponse(BaseModel):
    equation: str
    coefficients: list[int] = Field(description="全部系数（先反应物后生成物）")
    reactant_coefficients: list[int]
    product_coefficients: list[int]
    precision: int = Field(description="实际采用的小数位数")
    extent: str = Field(description="理论反应进度 ξ（mol，定精度小数字符串）")
    limiting_reactants: list[str] = Field(
        description="限制反应物（反应进度最小，可并列多个）"
    )
    products: list[ProductYield] = Field(description="各生成物理论产量")
    reactants: list[ReactantRemainder] = Field(description="各反应物投料/消耗/剩余")


class FormulaInferRequest(BaseModel):
    elements: dict[str, Decimal] = Field(
        description="各元素质量百分比（%，均须为正；总和与 100 的偏差不超过 tolerance）",
        examples=[{"H": "5.93", "C": "35.36", "O": "58.71"}],
    )
    tolerance: Decimal = Field(
        default=Decimal("0.1"),
        gt=0,
        description=(
            "允许的最大绝对偏差（质量百分点），同时用于百分比总和校验与"
            "摩尔质量相对偏差（%）；必须为正数，默认 0.1"
        ),
    )
    molar_mass: Decimal | None = Field(
        default=None,
        description=(
            "实测摩尔质量（g/mol）；提供时取经验式摩尔质量最接近的正整数倍"
            "确定分子式，相对偏差须不超过 tolerance%"
        ),
    )


class FormulaElementRow(BaseModel):
    element: str
    subscript: int = Field(description="经验式中该元素的下标（互质正整数）")
    given_percent: str = Field(description="请求给出的质量百分比（%，定精度小数字符串）")
    theoretical_percent: str = Field(
        description="由候选经验式反算的理论质量百分比（%，定精度小数字符串）"
    )
    deviation: str = Field(
        description="逐元素偏差：理论百分比 - 实测百分比（百分点，带符号）"
    )


class FormulaInferResponse(BaseModel):
    empirical_formula: str = Field(description="经验式（Hill 顺序）")
    subscripts: dict[str, int] = Field(description="各元素经验式下标（互质正整数）")
    formula_scale: int = Field(description="入选候选对应的摩尔比放大倍数（1～12）")
    total_atoms: int = Field(description="经验式总原子数")
    elements: list[FormulaElementRow] = Field(
        description="逐元素下标、理论百分比与偏差（按请求元素顺序）"
    )
    max_abs_deviation: str = Field(
        description="各元素绝对偏差的最大值（百分点，定精度小数字符串）"
    )
    empirical_molar_mass: str = Field(
        description="经验式摩尔质量（g/mol，定精度小数字符串）"
    )
    measured_molar_mass: str | None = Field(
        default=None, description="实测摩尔质量（g/mol），未提供时为 null"
    )
    molecular_multiplier: int | None = Field(
        default=None,
        description="分子式相对经验式的正整数倍（n），未提供摩尔质量时为 null",
    )
    molecular_formula: str | None = Field(
        default=None, description="分子式（Hill 顺序），未提供摩尔质量时为 null"
    )
    molecular_molar_mass: str | None = Field(
        default=None,
        description="分子式理论摩尔质量（g/mol），未提供摩尔质量时为 null",
    )
    molar_mass_relative_deviation_percent: str | None = Field(
        default=None,
        description="分子式摩尔质量与实测值的相对偏差（%），未提供摩尔质量时为 null",
    )
