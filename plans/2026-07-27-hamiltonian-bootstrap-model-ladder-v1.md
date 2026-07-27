# Hamiltonian Bootstrap Model Ladder

**文件状态：** In Progress  
**计划日期：** 2026-07-27  
**目标仓库：** `/home/hzxiaxz/Quantum-harness/Hamiltonian-Bootstrap`  
**实施路线（固定顺序）：** MG/显式二聚化 chain → cluster model 加场 → toric code 加场 → Kitaev honeycomb → Shastry–Sutherland → triangular J1-J2

## Objective

建立一条可重复、可诊断、按模型逐级放行的 Hamiltonian-bootstrap 研究路径，回答统一问题：**低阶局域 certificate 从 frustration-free/可解点出发，在扰动下如何失效；在固定 SDP budget 下，operator-adapted basis 能否比统一局域 basis 延长证书或数值界的有效区间。**

第一里程碑不是泛化或重写现有 optimized `GSB`，而是在当前 `Vector{UInt16}` Pauli word 之上建立最小共同的显式编译路径：调用者显式给出 Hamiltonian polynomial、basis sectors/seeds、symmetry generators 及其用途、linear state-optimality tests、PSD state-optimality basis、RDM regions、observables 与 certificate scope；程序只机械完成 Pauli products、`B†B`、`[H,A]`、论文式 PSD state-optimality entries、canonical moments、用户指定的 symmetry quotient/Fourier blocks，以及 JuMP constraints。编译和求解必须分离，使结构测试、budget 比较和数值诊断不依赖立即求解。

预期结果包括：

- 一个可检查的显式 relaxation specification 与 compiled artifact；
- 对 Ising 与 Heisenberg 旧入口保持回归兼容；
- 六个模型按固定顺序通过独立 stage gate，而不是一次性宣称通用；
- 每个模型都有 uniform-local 与 operator-adapted 两组固定-budget 实验；
- 数值求解结果、有效 strengthening 和严格 certificate 被明确区分；
- 不引入 RG，也不把 coarse graining 纳入本路线。

## Non-goals

- 不实现 auto symmetry detection，不从 Hamiltonian 自动猜测对称群。
- 不实现自动 basis selection、机器学习选基或自动补全 basis；basis enrichment 由实验配置显式声明。
- 不设计通用 qudit、fermion 或任意非交换代数 IR；第一阶段继续使用当前 `Vector{UInt16}` Pauli word 编码。
- 不把 PauliStrings.jl 设为第一阶段 blocker；它只可在核心路径稳定后作为可选输入 adapter，adapter 输出仍须落到本仓库 canonical Pauli polynomial。
- 不设计 serialization、schema migration、结果数据库或配置文件格式；实验 metadata 保持内存结构及终端/测试输出。
- 不重写现有 optimized `GSB`、Fourier 优化、RDM 8/9/10 或认证后端；它们保留为论文复现和 Heisenberg 回归后端。
- 不承诺为每个数值 SDP 解提供严格浮点后验认证；只有通过明确 certificate scope 和独立验证链的结果才称为严格证书。
- 不以“能量界较紧”替代 observable-specific basis adequacy；远距离 observable 必须单独评估。
- 不引入 RG、DMRG coarse-graining map 或多尺度重整化流程。

## Authoritative Design Decisions

### 1. 输入完全显式

调用者必须显式提供下列对象；编译器不得推断或猜测：

- Hamiltonian polynomial：canonical Pauli words、实/复系数、有限晶格展开方式和总能量/能量密度归一化；
- basis specification：具名 sectors、每个 sector 的 seeds、允许的显式 translation orbit 或用户列出的完整 words；
- symmetry generators：每个 generator 的作用，以及用途标签，例如 `moment_zero`、`moment_equality`、`basis_block`、`fourier_orbit` 或 `redundant_block_equivalence`；
- linear state-optimality test operators `A`；
- PSD state-optimality basis `\widetilde B`；
- RDM regions，包含有序 site 列表及是否允许用户声明的 block decomposition；
- observable polynomials；
- certificate scope，例如 `numerical_relaxation`、`solver_bound` 或 `rigorously_postvalidated`。

现有 `PauliSymmetryModel` 只有轴置换、全局 sign generators、translation/reflection 布尔值，不能承担多子格、方向相关轴变换或用途语义；其当前字段是明确的兼容边界，而非最终自动 symmetry API。`src/basic_function.jl:208-222`

### 2. 编译器只做机械生成

编译路径仅执行：

1. canonical Pauli multiplication，保留 `±1, ±im` 相位；
2. 从显式 basis 生成全部 `B_i†B_j` moment support；
3. 从显式 linear tests 生成 `H*A-A*H`；
4. 从显式 PSD basis 生成论文式条目
   \[
   K_{v,w}=\ell\!\left(vHw^\dagger-\tfrac12(Hvw^\dagger+vw^\dagger H)\right);
   \]
5. 为显式 RDM region 生成 Pauli expansion 所需 moments；
6. canonical moment 去重；
7. 应用用户指定且用途明确的 symmetry quotient/equality/block/Fourier 规则；
8. 生成 JuMP 变量、线性等式和 PSD constraints。

论文给出 linear 与 PSD state-optimality 的有限测试形式，并指出同一 symmetry 分块可用于 PSD state-optimality。`MinerU_markdown_2604.01555v1_2081746923717435392.md:590-616`

### 3. 四类结论不得混称

| 类别 | 含义 | 可声明内容 |
|---|---|---|
| 启发式 basis 选择 | 人工选择 sectors/seeds，以期在固定 budget 下改善界 | 只能称 basis policy 或实验选择；不改变约束的数学含义 |
| 严格等价约化 | 有证明的 symmetry quotient、moment equality、零 moment 或 Fourier block diagonalization | 与未约化 relaxation 等价；必须记录 generator、用途和适用 Hamiltonian |
| 有效 strengthening | RDM PSD、linear/PSD state-optimality、额外有效不等式 | 缩小可行域，但必须检查所需 moment support 完整；不能静默丢项 |
| 严格证书 | 经明确证书 scope 和独立后验验证得到的可审计结论 | 仅在认证流程成功后使用“严格/rigorous”；普通 solver objective 只报告数值下界候选或 solver bound |

论文的某些 Hamiltonian-only sign reductions 不是普通 Pauli 代数 automorphism，不能被泛化成自动群商；其正确性依赖论文中的 parity/real-block 论证。`MinerU_markdown_2604.01555v1_2081746923717435392.md:260-392`

### 4. 最小改动与文件策略

- 第一里程碑优先演化 `src/basic_function.jl` 与 `src/ising.jl`，并在 `test/runtests.jl` 增加结构与回归测试；不先拆分现有大文件。
- 通用 Pauli polynomial、basis specification、compiled artifact 和诊断类型优先放入现有 `src/basic_function.jl`，因为 Pauli 编码、reduction、symmetry 和 solver 参数已经在那里。现有 Pauli reduction 主入口位于 `src/basic_function.jl:323-336`。
- `src/ising.jl` 先改为新编译路径的薄模型 wrapper，同时保持 `ising_ground_state_bound` 与 `IsingMomentResult` 行为；当前 Ising 原型已经把 basis 与 moment-map 生成分成独立步骤。`src/ising.jl:11-56`
- `src/bound_gsp.jl` 不在第一阶段重构；它继续作为 optimized Heisenberg/论文复现后端。其当前函数把 basis、Fourier、约束、求解与结果提取放在同一次调用内。`src/bound_gsp.jl:1-2` `src/bound_gsp.jl:226-366` `src/bound_gsp.jl:548-732`
- 只有当第二个真实模型使用点已经出现、接口已被两处代码验证且共享代码确实必要时，才允许新建共享源码文件（候选名 `src/pauli_relaxation.jl`）；若现有文件仍清晰，则继续保留在现有文件。
- 模型配置先放在 `examples/ground_state.jl` 的新增具名函数或一个必要的模型阶梯 example 中；只有模型构造需要成为公共 API 时，才创建模型专属源码文件。现有 ground-state examples 已包含 chain 与 square J1-J2 调用模式。`examples/ground_state.jl:3-30` `examples/ground_state.jl:48-53`
- `src/QMBCertify.jl` 只在新公共类型/API 稳定后追加 include/export，并保持旧 export；当前入口和 include/export 全貌为 `src/QMBCertify.jl:1-43`。
- 测试先继续集中在 `test/runtests.jl`；仅当文件显著失去可读性时再拆分。当前全部覆盖只有两个 testset。`test/runtests.jl:1-34`

### 5. 数据与数值约定

- 第一阶段 canonical word 固定为 `Vector{UInt16}`，空 vector 是 identity；现有 Ising basis 使用同一约定。`src/ising.jl:11-33`
- 轴编码固定为 `1=X, 2=Y, 3=Z`，编码为 `3(site-1)+axis`；各向异性模型启用前必须增加测试，避免 certification helper 中 `z,x,y` 变量排列造成误用。核心 Ising 目标以 `UInt16[1]` 和 `UInt16[3,6]` 表示 `X₁` 与 `Z₁Z₂`。`src/ising.jl:108-119` certification 变量创建采用不同文本顺序。`src/certification/helpers.jl:67-73`
- Hamiltonian 容器名称可调整，但语义必须是显式 Pauli polynomial，而不是依赖 SU(2) 代表项的隐式展开。
- canonicalization 顺序必须确定；所有 dictionary 输出在形成 compiled artifact 前按稳定 key 排序，不依赖 hash iteration。
- 不再沿用 LSO 的随机投影去重作为新路径的一部分；现有随机选择发生在 `filter_mons`。`src/bound_gsp.jl:818-872`
- 新路径所有 state-optimality 系数至少使用 `Float64`；不得复制现有 PSO `Float16` 中间系数。现有降精度位置是 `src/bound_gsp.jl:756-758`。

## Current-state Map

### Project Structure Summary

| 区域 | 当前现实 | 实施含义 |
|---|---|---|
| 包入口 | JuMP/Mosek 等依赖、exports、`qmb_data` 与主体 includes 集中在模块文件。`src/QMBCertify.jl:1-43` | 只在公共接口稳定后修改；必须保持 `GSB`、Ising 和 certification exports |
| Heisenberg basis | `get_basis` 是 chain/square、label、degree、extra 的手工模板生成器。`src/basic_function.jl:10-123` | 作为旧后端兼容合同，不当作任意模型 basis generator |
| Pauli algebra/symmetry | reduction、`PauliSymmetryModel`、chain/square canonicalization 已存在。`src/basic_function.jl:143-348` | 复用 Pauli 乘法；新显式 symmetry 用途层不得破坏旧默认行为 |
| Ising prototype | `_pauli_basis` 和 `_ising_moment_data` 已形成 basis → moment map 原型。`src/ising.jl:11-56` | 第一里程碑从这里抽象出显式 polynomial/basis/compile 路径 |
| Ising solve | 目标 Hamiltonian 和 Mosek 求解仍在一个函数中。`src/ising.jl:71-132` | 需要 compile/solve 分离，并用 wrapper 保持旧 API |
| Fourier/GSB | chain/square Fourier、Gram blocks、LSO/PSO、RDM、目标、求解和提取集中实现。`src/bound_gsp.jl:38-221` `src/bound_gsp.jl:226-599` | 不先重写；只提取行为测试和后续经验证的公式 |
| RDM | 只有硬编码 8/9/10-site Heisenberg blocks。`src/rdm_positivity.jl:1-50` | 新显式 RDM region 先走小 region 的机械 Pauli expansion，不泛化旧硬编码后端 |
| 测试 | 仅 symmetry 与 Ising 两个 exact endpoints。`test/runtests.jl:1-34` | 第一 gate 必须先补纯编译、确定性与新旧回归 |
| 论文 sparse basis | 一维连续 strings 与可选远距 pairs、二维局域模板有明确描述。`MinerU_markdown_2604.01555v1_2081746923717435392.md:189-201` | 作为显式 basis policy 示例，不作为自动选基规则 |
| 论文 symmetry | sign、real、translation/Fourier、axis permutation、reflection/D4 有严格推导。`MinerU_markdown_2604.01555v1_2081746923717435392.md:212-560` | 每项约化必须由用户显式启用并标注用途/适用条件 |
| 论文 strengthening | RDM、linear/PSD state-optimality 给出公式。`MinerU_markdown_2604.01555v1_2081746923717435392.md:576-616` | 新编译器机械生成并做 closure 检查 |
| 论文 benchmark | J1-J2、MG energy/observables 有数值基准。`MinerU_markdown_2604.01555v1_2081746923717435392.md:645-724` | 作为第一模型的物理验收，不把能量闭合等同于 observable 闭合 |

### Table 2 的权威解释

Table 2 中的 **31** 是 **`N=100, d=4, r=1` 在利用论文全部所列结构约化后，最大的单个 PSD block 的矩阵阶数**；它不是约束数、moment 数、PSD block 数、变量数或完整 moment matrix 的尺寸。论文明确把表描述为 maximal block sizes；偶数 `d` 的公式为 `(3^(d+1)+5)/8`，代入 `d=4` 得 31。`MinerU_markdown_2604.01555v1_2081746923717435392.md:562-574`

同一表的结构回归链应记录：原始 8,127,090,301；Pauli equalities 后 322,029,976；locality sparsity 后 12,001；全部 symmetry 后最大块 31。这里 12,001 是 sparse basis size，而 31 是最大块尺寸，二者不得混用。`MinerU_markdown_2604.01555v1_2081746923717435392.md:562-574`

## Prioritized Challenges and Risks

1. **P0—严格性标签错误。** 普通 solver objective 若被称为“严格 certificate”，会使整个研究结论不可审计。优先在数据类型和输出中固定 certificate scope。
2. **P0—缺失 moment 被静默忽略。** `B†B`、commutator、PSD state-optimality 与 RDM 都可能产生 basis 外 support；必须编译时报错或按显式 scope 分类，不能补零。
3. **P0—错误 symmetry quotient。** 带场、bond-dependent 或 patterned Hamiltonian 不具备现有 Heisenberg symmetry；任何自动复用都可能改变问题。故 symmetry 完全显式且用途化。
4. **P1—编译非确定性与不可比较 budget。** 当前 LSO 随机投影和 hash 顺序会破坏结构回归；先建立稳定排序、exact row signature 与 diagnostics。
5. **P1—编译/求解耦合。** 当前 Ising 与 `GSB` 都会在构造后立即求解。`src/ising.jl:79-124` `src/bound_gsp.jl:580-599` 先分离才能做无许可证结构测试和跨模型 budget 控制。
6. **P1—数值病态。** 论文报告 state-optimality 在一维部分区间和二维模型产生数值问题。`MinerU_markdown_2604.01555v1_2081746923717435392.md:783-784` 必须逐层启用并保留降级路径。
7. **P2—几何多样性。** 二聚化、多子格、边 qubit、honeycomb 轴-键绑定、orthogonal dimers 和 triangular primitive vectors 不能塞进当前 chain/square 布尔分支。
8. **P2—basis 与 observable 失配。** 论文二维最远距离关联随尺寸显著恶化，明确归因于高度局域 basis。`MinerU_markdown_2604.01555v1_2081746923717435392.md:746-756` 每个模型必须有 observable-adapted ablation。
9. **P3—过早抽象。** 预建通用 lattice/qudit/serialization 层会延迟可验证结果；共享源码只在第二真实使用点出现后抽取。

## Implementation Plan

### Phase 0 — Baseline Lock and Terminology

- [x] **0.1 [Done] 固定旧行为清单与测试基线。** 修改位置：`test/runtests.jl:1-34`。记录现有 Ising 两个 exact endpoints、现有 symmetry canonical representatives、`IsingMomentResult` 字段和 `GSB` 公开调用/返回形状；对不应在无许可证环境运行的 solver tests 使用现有测试策略而不改变数学期望。理由：后续重构必须能区分预期 API 演化与回归。
- [x] **0.2 [Done] 增加 Pauli 编码/乘法的纯结构测试。** 修改位置：`test/runtests.jl:1-34`；实现位置仍为 `src/basic_function.jl:143-206`。覆盖空 word、`XX=I`、`XY=iZ`、`YX=-iZ`、异站点排序、`realify` 差异和 `UInt16` 最大 site 检查。理由：后续六个各向异性/多体模型都依赖这一层，轴语义错误必须在进入模型前暴露。
- [x] **0.3 [Done] 固定术语与结果 scope。** 修改位置：优先 `src/basic_function.jl` 中新增最小 enum/symbol validation；公共后再修改 `src/QMBCertify.jl:15-20`。保证 diagnostics 和 solve result 显式区分 heuristic basis、equivalent reduction、strengthening、numerical solver bound、rigorously postvalidated certificate。理由：防止研究输出越权声明。

**Phase 0 完成定义：** 纯 Pauli 测试确定性通过；旧 symmetry 与 Ising tests 不变；任何结果都不能在 scope 缺失时标记为 strict certificate。

### Phase 1 — Minimal Explicit Compile Path (First Milestone)

- [ ] **1.1 [Not Started] 定义显式 Pauli polynomial/Hamiltonian term 容器。** 修改位置：优先 `src/basic_function.jl`，仅在第二使用点证明确需拆分时迁至新 `src/pauli_relaxation.jl`；稳定后更新 `src/QMBCertify.jl:15-43`。容器至少保存 canonical `Vector{UInt16}` word 与 coefficient，构造时合并重复项、删除精确零项、稳定排序、验证 Hermiticity/允许的复系数及 site 编码范围；Hamiltonian 与 observable 使用同一 polynomial 值对象，但在 specification 中角色分离。理由：替代 `supp`/`coe` 并行数组的隐式语义，同时不设计大型 IR。
- [ ] **1.2 [Not Started] 定义显式 basis specification。** 修改位置同 1.1。表示具名 sectors 与显式 seeds/orbits；每个 sector 记录 `name`、ordered words、PSD 角色、可选用户声明的 translation orbit metadata。禁止从 Hamiltonian 自动生成或猜测 sectors。理由：使 uniform-local 与 operator-adapted basis 可在相同 budget 下精确重放。
- [ ] **1.3 [Not Started] 定义 relaxation specification。** 修改位置同 1.1。聚合 Hamiltonian、basis sectors、用途化 symmetry declarations、linear tests、PSD state-optimality basis、RDM regions、observables、normalization 与 certificate scope；构造时验证所有必需字段。理由：清空上下文后，单个对象即可完整描述要编译的问题。
- [ ] **1.4 [Not Started] 从任意显式 basis 自动生成 moment support。** 修改位置：演化 `src/ising.jl:35-56` 的 `_ising_moment_data`，模型无关逻辑稳定并有第二使用点后放入 1.1 所在共享边界。对每个 sector 机械生成全部 `reverse(B_i); B_j`，调用 canonical Pauli multiplication，保留相位、Hermitian 对应关系、零 moment 原因，并生成稳定的 canonical moment index。理由：这是任意模型共同的最小 moment compiler。
- [ ] **1.5 [Not Started] 编译 linear state optimality。** 修改位置：通用编译边界；旧 `GSB` 逻辑 `src/bound_gsp.jl:370-415` 仅作行为参考，不直接重写。对用户显式给出的每个 `A` 生成完整 `[H,A]`，canonicalize、合并系数、形成确定性线性行；使用 exact canonical row signature 去零和去重，不用随机投影。理由：满足通用多体/各向异性 Hamiltonian，而非现有二体 SU(2) 假设。
- [ ] **1.6 [Not Started] 编译论文式 PSD state-optimality entries。** 修改位置：通用编译边界；参考现有 `PSDstate_entry` `src/bound_gsp.jl:756-786`，但使用 `Float64` 或更高精度输入，不复制 `Float16`。仅从用户显式 `\widetilde B` 生成完整 `K_{v,w}`，保证 Hermitian 配对；任何缺失 canonical moment 均进入 deterministic support error。理由：把 strengthening 与 basis 选择解耦，并忠实实现论文公式。`MinerU_markdown_2604.01555v1_2081746923717435392.md:598-616`
- [ ] **1.7 [Not Started] 编译显式小区域 RDM。** 修改位置：优先通用编译边界；旧固定后端保留在 `src/rdm_positivity.jl:1-50`。按用户 region 的 `[I,X,Y,Z]^k` Pauli expansion 机械生成所需 moments，默认不自动猜测 magnetization blocks；只有用户显式声明并经测试的 block decomposition 才应用。理由：toric star/plaquette、honeycomb hexagon和三角 cluster 需要几何显式 region。
- [ ] **1.8 [Not Started] 应用用户指定的 symmetry 用途。** 修改位置：扩展 `src/basic_function.jl:208-348` 周围的显式声明层，但保持 `PauliSymmetryModel` 与 `reduce!` 旧方法。分别实现 moment-zero/equality、basis-block、用户给定 orbit canonicalization、显式 Fourier orbit；每项检查 Hamiltonian invariance 由调用者声明且编译器做机械一致性校验，失败则拒绝编译。不得自动发现 generator。理由：把严格等价约化与经验 basis policy 分开。
- [ ] **1.9 [Not Started] 建立 compiled artifact 与 compilation diagnostics。** 修改位置：通用编译边界。artifact 至少包含稳定 moment table、每个 PSD block 的 affine entries、linear rows、objective rows、constraint provenance 和 source specification reference；diagnostics 至少包含 raw/canonical basis size、raw `B†B` entries、zeroed/equated moments、scalar moment count、linear row raw/zero/duplicate/independent counts、PSD block count/dimensions、最大 block size、real-embedded dimensions、RDM/state-optimality 增量、missing-support 明细和 deterministic fingerprint。理由：无需求解即可验收结构与 fixed-budget 公平性。
- [ ] **1.10 [Not Started] 分离 compile 与 solve。** 修改位置：`src/ising.jl:71-132` 及通用边界。`compile_relaxation(spec)` 不创建绑定特定 solver 的最终结果且绝不调用 `optimize!`；`build_jump_model(compiled; optimizer...)` 只机械创建 JuMP constraints；`solve_relaxation(...)` 负责 optimizer、状态和数值诊断。旧 `ising_ground_state_bound` 变为 compile+solve wrapper，保持参数、返回结果和端点数值。理由：支持纯结构 CI、solver injection 和重复求解 observables。
- [ ] **1.11 [Not Started] 加入 deterministic support/constraint checks。** 修改位置：通用编译边界和 `test/runtests.jl`。覆盖 Hamiltonian/objective moment closure、`B†B` Hermiticity、linear row residual structure、PSD state-optimality Hermiticity、RDM trace/identity support、symmetry orbit coefficient consistency、重复 constraints、空 PSD block、非有限系数和稳定 fingerprint。所有失败提供来源（basis pair、test operator、Hamiltonian term、region）。理由：禁止静默补零或漏项。
- [ ] **1.12 [Not Started] 完成新旧 Ising/Heisenberg 回归。** 修改位置：`src/ising.jl`、`test/runtests.jl`，必要时只为测试调用现有 `GSB`。Ising 新 wrapper 复现 `J=0` 与 `h=0` endpoints；Heisenberg 使用现有 `get_basis` 与 `GSB` 小规模配置，核对旧 basis 轨道布局、support、目标值/状态和 `qmb_data` 返回形状。不得重写 optimized GSB。理由：第一共同路径上线前必须证明未破坏论文复现后端。
- [ ] **1.13 [Not Started] 复现 Table 2 结构 benchmark。** 修改位置：`test/runtests.jl` 或一个明确的非默认重型结构测试。显式构造论文一维 `N=100,d=4,r=1` sparse basis policy，验证 sparse size 12,001 及全部用户指定结构约化后的最大 PSD block size 31；测试名称必须写明“31 is max PSD block dimension, not constraint count”。理由：验证 basis、sector、translation/Fourier 和 diagnostics 的共同语义。`MinerU_markdown_2604.01555v1_2081746923717435392.md:562-574`

**Phase 1 完成定义：** 可在不求解的情况下从完全显式 specification 得到确定性 compiled artifact；任意 basis 的 `B†B`、通用 linear/PSD state-optimality 与小 RDM region 均有 closure 检查；Ising 旧入口与 Heisenberg `GSB` 回归通过；Table 2 的 12,001 与最大块 31 被正确区分并复现。

### Phase 2 — MG and Explicitly Dimerized Chain

- [ ] **2.1 [Not Started] 建立模型显式 Hamiltonian。** 修改位置：优先 `examples/ground_state.jl` 的具名构造与测试 helper；只有成为公共 API 才新增模型专属源码。采用 Pauli 约定
  \[
  H=\tfrac14\sum_i J_1[1+(-1)^i\delta](X_iX_{i+1}+Y_iY_{i+1}+Z_iZ_{i+1})+\tfrac{J_2}{4}\sum_i(X_iX_{i+2}+Y_iY_{i+2}+Z_iZ_{i+2}).
  \]
  显式展开偶/奇 bond classes，不使用 SU(2) 代表项隐式乘 3；`δ=0,J2/J1=1/2` 是 MG anchor，`δ=1,J2=0` 是解耦 dimer anchor。理由：第一个模型同时验证显式 polynomial 与 period-2 geometry。
- [ ] **2.2 [Not Started] 声明两组 basis policy。** uniform-local 使用显式连续 Pauli strings（先从长度 1–2，再按 budget 加至 3–4）；operator-adapted 在相同 budget 下优先加入 strong-dimer singlet-channel 等价 Pauli words、相邻两个 dimer 的 products、J2 跨 dimer pairs 和目标 observable support。sectors 与 seeds 全部列出，不自动生成猜测。理由：直接检验“从 frustration-free/可解点扰动时 operator-adapted basis 是否更耐用”。
- [ ] **2.3 [Not Started] 声明 symmetry 和 strengthening。** `δ=0` 可显式使用一站点 translation；`δ≠0` 只允许两站点 translation；reflection、全局 sign 和轴 permutation 必须按参数点显式启用并验证 Hamiltonian invariance。linear tests 与 PSD basis 从小局域 dimer operators 显式列出，逐层启用。理由：防止当前 `translation::Bool` 把 period-2 问题错误商掉。
- [ ] **2.4 [Not Started] 增加 observables 与 benchmark。** observables 至少含能量密度、强/弱 bond energy、dimer order、`C(1)=<X_iX_{i+1}>/4`、`C(2)=<X_iX_{i+2}>/4`。MG benchmark：`e0=-3/8`；`C(1)` interval 覆盖 `-1/8`，`C(2)` interval 覆盖 0。论文基准与 interval 位于 `MinerU_markdown_2604.01555v1_2081746923717435392.md:645-699`。
- [ ] **2.5 [Not Started] 执行扰动扫描与停止规则。** 从 `(δ=0,J2/J1=0.5)` 沿 `J2/J1` 双向扰动，并从 `(δ=1,J2=0)` 降低 `δ`；每点比较相同 max-block/scalar-moment budget。若 MG 能量不能闭合、解耦 dimer 端点不精确、period-2 symmetry 检查失败，或 adapted basis 在至少一个非零扰动点未优于 baseline 且诊断无法解释，则停在本 stage，不进入 cluster。

**Phase 2 完成定义：** 两个 anchor 通过；uniform/adapted fixed-budget 曲线、observable intervals、constraint provenance 和失效阈值可重放；MG 结果与旧 `GSB` 论文后端交叉核对。

### Phase 3 — Cluster Model with Field

- [ ] **3.1 [Not Started] 建立显式 cluster Hamiltonian。** 预期模型专属位置：先在 example/test helper；若公共使用，新增 `src/cluster_chain.jl` 并更新 `src/QMBCertify.jl:35-43`。采用明确约定
  \[
  H=-J\sum_i Z_{i-1}X_iZ_{i+1}-h_x\sum_iX_i-h_z\sum_iZ_i,
  \]
  首轮固定 `h_z=0`，另以 `h_z` 作为显式 symmetry-breaking stress test。理由：验证三体 Hamiltonian、场项和非 Heisenberg symmetry。
- [ ] **3.2 [Not Started] 声明 basis/observables。** baseline 为连续局域 strings；adapted seeds 包含 cluster stabilizers `Z-X-Z`、相邻 stabilizer products、field operators、稳定子端点 string fragments。observables 为能量、`<Z_{i-1}X_iZ_{i+1}>`、`<X_i>`、有限长度 cluster string order；目标 string 的端点和路径必须显式进入 adapted basis。
- [ ] **3.3 [Not Started] 声明 symmetry/strengthening。** 不复用 Heisenberg 轴置换；由用户显式列出实际全局/子格 symmetry generators及其用途。先 baseline PSD，再 linear commutators，再小 PSD state-optimality basis；每层记录 conditioning。理由：测试用途化 symmetry 和通用三体 commutator。
- [ ] **3.4 [Not Started] Benchmark 与停止规则。** `h_x=h_z=0` commuting stabilizer 点必须达到精确能量并固定 stabilizer expectation；小系统用 exact diagonalization 对照 field sweep。若三体项 closure、string observable support 或 symmetry-breaking `h_z` 的 generator rejection 任一失败，停止；若 adapted basis 在固定 budget 下不能延长 stabilizer/string interval 的有效区间，先完成 basis ablation 后再决定是否进入 toric。

**Phase 3 完成定义：** commuting anchor、场扰动、小系统 ED、string observable 和 symmetry-breaking negative test 全部通过；这是通用 polynomial/compiler 的第二真实模型用点。此时才评估是否把共享编译代码从现有文件移入 `src/pauli_relaxation.jl`。

### Phase 4 — Toric Code with Field

- [ ] **4.1 [Not Started] 建立 edge-qubit 周期几何和 Hamiltonian。** 预期模型专属文件 `src/toric_code.jl` 仅在模型进入公共 API 时创建；初版明确 `(cell_x,cell_y,edge_orientation)→site`，生成 star `A_s=∏X_e`、plaquette `B_p=∏Z_e` 与
  \[
  H=-J_s\sum_sA_s-J_p\sum_pB_p-h_x\sum_eX_e-h_z\sum_eZ_e.
  \]
  不把当前每 cell 一个 site 的 square 分支复用为 toric geometry。理由：验证多 qubit unit cell 与四体 stabilizers。
- [ ] **4.2 [Not Started] 声明 basis/regions/observables。** baseline 为半径受限局域 words；adapted seeds 包含 stars、plaquettes、相邻 stabilizer products、field-dressed stabilizers，以及显式短 Wilson/'t Hooft loop segments。RDM regions 显式列出单 star、单 plaquette 和相邻 star-plaquette union。observables 为能量、`<A_s>`、`<B_p>`、短闭环与开放 string endpoints。
- [ ] **4.3 [Not Started] 声明 symmetry 与 sector scope。** translation、edge-orientation exchange、point-group 或 stabilizer-related equalities均由用户明确列出；local stabilizers不能在加场后继续当作无条件 symmetry quotient。拓扑 sector 若未由有限局域 constraints 固定，必须在 certificate scope 中声明“不区分全局 sector”。理由：避免把可解点结构错误外推到扰动区。
- [ ] **4.4 [Not Started] Benchmark 与停止规则。** `h_x=h_z=0` commuting-projector 点必须精确；纯单方向场先做小 torus ED。若 edge indexing、star/plaquette commuting tests、loop support、RDM region 或 sector scope 不完整，停止；若数值 state-optimality 导致不稳定，允许退回 baseline/linear-only，但必须记录降级，不能删除失败点。

**Phase 4 完成定义：** commuting toric anchor、小场 ED、局域 stabilizer 与短 loop intervals 通过；compiled diagnostics 能区分局域块与全局 sector 范围。

### Phase 5 — Kitaev Honeycomb

- [ ] **5.1 [Not Started] 建立 honeycomb 两子格几何和 bond-dependent Hamiltonian。** 预期 `src/kitaev_honeycomb.jl`；采用
  \[
  H=-J_x\sum_{\langle ij\rangle_x}X_iX_j-J_y\sum_{\langle ij\rangle_y}Y_iY_j-J_z\sum_{\langle ij\rangle_z}Z_iZ_j-\sum_i(h_xX_i+h_yY_i+h_zZ_i).
  \]
  显式列出 bond types 与轴绑定。若 toric 与 Kitaev 已证明共享相同 cell/site/wrap 协议，才抽最小 `src/unit_cell_lattices.jl`；否则保持模型专属。理由：验证空间操作与 Pauli 轴联动，而非纯内部轴 symmetry。
- [ ] **5.2 [Not Started] 声明 basis/observables。** baseline 为局域 bond words；adapted seeds 包含 x/y/z bond terms、hexagon flux operator、相邻 bond products、field insertions和 flux-neighborhood clusters。observables 为能量、各 bond energy、plaquette flux、短 spin correlations。
- [ ] **5.3 [Not Started] 声明 symmetry 与可解 benchmark。** 零场使用用户明确给出的 translation/point-group+axis 联合作用；不得把 `S3` 轴置换独立应用于各向异性 couplings。`h=0` 与已知精确有限-size/自由费米子能量或 ED 对照；场扰动从小场开始。
- [ ] **5.4 [Not Started] 停止规则。** 若联合空间-轴 generator 不能通过逐 term invariance，若 flux operator canonicalization 错误，或零场 benchmark 超出预设数值容差，停止；只有 adapted flux basis 在固定 budget 下至少改善 flux/能量有效区间之一，且结果不依赖错误 symmetry，才进入下一模型。

**Phase 5 完成定义：** 零场可解 benchmark、bond/flux observables、联合 symmetry negative/positive tests 和小场扫描可重放。

### Phase 6 — Shastry–Sutherland

- [ ] **6.1 [Not Started] 建立 patterned orthogonal-dimer Hamiltonian。** 预期 `src/shastry_sutherland.jl`；显式生成 square-edge `J` bonds 与选定正交对角 `J'` dimer bonds，使用扩大 unit cell，不把它当作 uniform square J1-J2。Pauli 展开每个 `S_i·S_j` 为 `(XX+YY+ZZ)/4`。理由：验证 patterned bond classes 和显式 dimer product anchor。
- [ ] **6.2 [Not Started] 声明 basis/observables。** baseline 为局域 square clusters；adapted seeds 为 exact dimer bond operators、相邻 orthogonal dimers、inter-dimer edge terms及目标 plaquette/dimer correlations。observables 为能量、dimer bond energy、inter-dimer correlation、dimer order和小 plaquette quantity。
- [ ] **6.3 [Not Started] Benchmark 与扰动。** 从 `J/J'=0` 解耦 dimer product 点开始增大 `J/J'`，用小 cluster ED 和已知 dimer-product 能量交叉核对。symmetry 只使用 patterned lattice 实际保留的 translation/point group，全部显式列出。
- [ ] **6.4 [Not Started] 停止规则。** 若扩大 unit-cell bond enumeration、dimer product anchor 或 symmetry orbit counts 不匹配，停止；若 fixed-budget adapted basis 不能比 uniform basis 延长 dimer observable 的窄 interval 区间，先完成 seeds/budget ablation，不能以增加无界 budget 通过 gate。

**Phase 6 完成定义：** 解耦 dimer anchor、至少一段非零 `J/J'` 扫描、ED 对照和 fixed-budget basis 对比完成。

### Phase 7 — Triangular J1-J2

- [ ] **7.1 [Not Started] 建立 triangular primitive-cell Hamiltonian。** 预期 `src/triangular_j1j2.jl`；显式枚举三组无向 J1 bond directions 与正确 J2 shell，使用 primitive coordinates 与周期 wrap，Pauli 展开 `S_i·S_j`。不复用 square `slabel`/D4 逻辑。理由：这是最终非二分挫折 stress test。
- [ ] **7.2 [Not Started] 声明现实 benchmark scope。** triangular antiferromagnetic J1-J2 不强行宣称存在与前述模型相同的 frustration-free anchor；以小周期 torus exact diagonalization 为权威起点，并将其明确标为“finite-size exactly benchmarked stress stage”，不是热力学严格可解点。统一研究问题在此检验从前五个 anchor 学到的显式 operator-adapted policy 是否仍有效，而不是伪造可解点。
- [ ] **7.3 [Not Started] 声明 basis/observables。** baseline 为固定半径局域 words；adapted seeds 包含三方向 J1/J2 bond energies、elementary triangle products、rhombus clusters和目标 structure-factor 位移 pairs。observables 为能量、分方向 bond correlations、triangle/rhombus quantities、有限尺寸 structure factors。
- [ ] **7.4 [Not Started] 声明 symmetry/strengthening。** translation 和 triangular point-group generators 显式给出；SU(2) 只在 Hamiltonian 参数确实保持时启用。按 baseline → linear → small PSD state-optimality → explicit RDM regions 分层，二维 state-optimality 数值失败是预期风险而非自动删除理由。论文明确报告二维 state-optimality 数值问题。`MinerU_markdown_2604.01555v1_2081746923717435392.md:783-784`
- [ ] **7.5 [Not Started] 停止规则。** bond shell/point-group orbit 与 ED Hamiltonian 必须逐项一致；若 baseline 已超过 budget、state-optimality 无法稳定或 observable support closure 失败，则保留已通过的较低层并停止扩张。只有在至少两个非平凡 `J2/J1` 点完成 fixed-budget baseline/adapted 对照后，模型阶梯才算完成。

**Phase 7 完成定义：** 小 torus ED、几何/symmetry 结构测试、至少两个耦合点的能量与 observable intervals、数值降级记录和 fixed-budget comparison 完成。

### Phase 8 — Optional Adapter and Certification Handoff

- [ ] **8.1 [Not Started] 评估 PauliStrings.jl 输入 adapter。** 仅当六模型核心路径稳定且存在真实外部输入需求时实现；adapter 只负责转成 canonical Pauli polynomial，并与原生 `Vector{UInt16}` 输入产生相同 deterministic fingerprint。理由：避免把外部包变成第一阶段 blocker。
- [ ] **8.2 [Not Started] 评估严格后验认证范围。** 复用现有 Gram rationalization/认证思路前，先明确 compiled artifact 中哪些块和 multipliers 足以重建严格证书；未覆盖的新 state-optimality/RDM/Fourier 类型保持 `solver_bound` scope。现有 certification 依赖 `qmb_data` 的 Gram/basis/support。`src/certification/energy_cert.jl:3-21`
- [ ] **8.3 [Not Started] 保持 optimized GSB 论文复现后端。** 用新 diagnostics 与旧 `GSB` 在 Heisenberg/MG benchmark 上并排报告，但不要求两个后端内部结构相同。理由：新通用路径优先可审计，旧路径优先规模和论文复现。

**Phase 8 完成定义：** 可选 adapter 不改变核心 IR；严格认证只覆盖明确支持的 compiled block 类型；旧 GSB 仍可独立运行。

## Model Specifications and Stage Gates

| Stage | Hamiltonian anchor → perturbation | Explicit basis comparison | Primary observables | Benchmark | Gate / stop condition |
|---|---|---|---|---|---|
| MG/dimer chain | MG `J2/J1=1/2,δ=0` 与 decoupled dimers `δ=1,J2=0` → `J2,δ` 扰动 | contiguous local vs dimer/J2 adapted | energy, strong/weak bonds, dimer order, C(1), C(2) | `e=-3/8`; C(1) covers `-1/8`; C(2) covers 0; old GSB | anchor、period-2 symmetry 或 fixed-budget comparison 失败即停 |
| Cluster+field | commuting cluster `h=0` → `h_x`，再以 `h_z` 破 symmetry | contiguous local vs stabilizer/string adapted | energy, stabilizer, field polarization, string order | exact commuting point + small ED | 三体 closure、symmetry-breaking rejection、string support 任一失败即停 |
| Toric+field | commuting `h_x=h_z=0` → one-axis then two-axis fields | local words vs star/plaquette/loop adapted | energy, star, plaquette, short loops | exact stabilizer point + small torus ED | edge geometry、commutation、sector scope、loop support 不完整即停 |
| Kitaev honeycomb | exact zero field → small field | local bonds vs flux/bond adapted | energy, x/y/z bonds, flux | exact/ED zero-field finite lattice | joint spatial-axis invariance 或 flux canonicalization 失败即停 |
| Shastry–Sutherland | decoupled orthogonal dimers `J/J'=0` → larger `J/J'` | square-local vs dimer/plaquette adapted | energy, dimer/inter-dimer correlations | exact dimer product + ED | patterned bonds、unit cell 或 dimer anchor 失败即停 |
| triangular J1-J2 | finite-size ED reference → selected `J2/J1` | local radius vs bond/triangle/rhombus/observable adapted | energy, directional correlations, structure factors | small torus ED; no false FF claim | geometry/orbits/ED mismatch，或 budget/conditioning 不允许至少两点比较即停 |

## Cross-model Experimental Matrix

每个允许求解的模型/参数点都必须运行下列显式配置；不允许因数值不利而从结果表中删除失败配置：

| Axis | Levels |
|---|---|
| Basis policy | `uniform_local`; `operator_adapted` |
| Budget tier | 以 `(max PSD block dimension, total PSD packed entries, scalar moment count, linear constraint count)` 的预声明上限定义 `small/medium`；比较时至少 max block 与 scalar moments 不得超过同 tier 上限 |
| Strengthening | baseline moment PSD; +linear SO; +small PSD SO; +explicit RDM；若上层失败保留下层 |
| Symmetry | none/minimal correctness baseline; declared equivalent quotient; declared quotient+Fourier（适用时） |
| Observable | energy；模型局域 order/stabilizer/bond；非局域或长程 observable |
| Perturbation | anchor；small；intermediate；first-failure neighborhood |
| Result scope | compile-only；numerical solver bound；rigorously postvalidated（仅支持时） |

固定-budget 主指标：

- energy lower-bound gap 或与 ED/known upper bound 的 bracket；
- observable interval width `UB-LB`；参考值为零时禁止使用相对 gap；
- 首个超过预声明有效阈值的扰动参数，作为“有效区间终点”；
- adapted 相对 uniform 的有效区间延长量；
- compilation diagnostics、wall-clock（仅作机器内参考）、solver iterations/status、residual 与最小 PSD eigenvalue；
- strengthening 的边际收益除以新增 scalar moments/PSD packed entries。

二维远距离 observable 必须单列结果，不能用局域能量替代。论文 square-lattice 最远距离关联从 `L=4` 到 `L=16` 的 interval 显著恶化，并把原因归于局域 basis 与远距离目标失配。`MinerU_markdown_2604.01555v1_2081746923717435392.md:746-756`

## Numerical Pathology Diagnostics

- [ ] **N.1 [Not Started] 编译级尺度诊断。** 对每个 affine row/block 报告 coefficient `min/max`、零行、重复行、近比例行候选、moment 使用次数和 block sparsity；对 objective 与 Hamiltonian 记录归一化尺度。
- [ ] **N.2 [Not Started] 求解级状态诊断。** 结构化记录 MOI termination/primal/dual status、objective/bound、primal/dual residual、relative gap、迭代数、solver warning；不可用字段明确为 `missing`，不伪造零值。
- [ ] **N.3 [Not Started] PSD 后验诊断。** 对返回的每个 Gram/state-optimality/RDM block 计算对称化后的最小 eigenvalue、最大 eigenvalue、估计 condition ratio 和违反量；实嵌入块同时核对复 Hermitian 重建一致性。
- [ ] **N.4 [Not Started] 分层隔离病因。** 按 baseline → linear SO → PSD SO → RDM 的顺序增量编译/求解；失败时保持同 basis/symmetry/budget，只移除最后新增层，定位是 support、scaling、redundancy 还是 cone conditioning。
- [ ] **N.5 [Not Started] 数值降级规则。** `ALMOST_OPTIMAL`、数值失败或残差超阈值不得自动当作通过；可保留结果用于诊断，但 scope 降为 `numerical_diagnostic`。禁止为得到好表格而静默移除 constraint blocks。
- [ ] **N.6 [Not Started] 精度规则。** 新编译路径不使用 `Float16`；对近零 canonical coefficients 使用显式 tolerance policy 并记录被删除项，默认优先 exact integer/Gaussian-integer Pauli phases与输入 Float64 系数分离。

论文已确认：一维 `0.1≤J2≤0.9` 的 PSD state-optimality 会有数值问题，二维 J1-J2 中 linear 与 PSD state-optimality 都可能出问题；因此失败是实验结果的一部分。`MinerU_markdown_2604.01555v1_2081746923717435392.md:783-784`

## Verification Criteria

- 编译同一 specification 两次得到完全相同的 ordered moments、constraints、block dimensions 和 deterministic fingerprint。
- Hamiltonian、observable、`B†B`、linear commutator、PSD state-optimality 与 RDM 所需 support 均被显式纳入或产生带 provenance 的编译错误；无静默缺项。
- `compile_relaxation` 不创建 solver side effect、不调用 `optimize!`；同一 compiled artifact 可用于多个 optimizer/observable solve。
- 旧 `ising_ground_state_bound` 两个 exact endpoints 与当前容差一致。`test/runtests.jl:18-34`
- 旧 Heisenberg basis 的轨道长度/顺序合同和至少一个小规模 `GSB` 结果保持；现有手工 basis 范围为 `src/basic_function.jl:10-123`。
- `N=100,d=4,r=1` 结构 benchmark 得到 sparse basis 12,001 和最大 PSD block size 31，并明确 31 不是 constraint count。`MinerU_markdown_2604.01555v1_2081746923717435392.md:562-574`
- 六个模型严格按固定路线通过 stage gate；未通过的模型不得被后续模型的局部成功掩盖。
- 每个模型至少有一个 anchor/ED benchmark、一个非零扰动点、一个局域 observable 和一个 operator-adapted fixed-budget comparison。
- 对 reference observable 为零的情况使用绝对 interval width；MG `C(2)=0` 是强制测试。`MinerU_markdown_2604.01555v1_2081746923717435392.md:689-699`
- 所有 symmetry reduction 都能追溯到用户声明 generator、用途和 Hamiltonian invariance check；没有 auto symmetry。
- 所有 basis words 都能追溯到用户声明 sector/seed/orbit；没有自动 basis 猜测。
- 数值失败点仍保留 diagnostics；constraint 降级是显式实验维度。
- 只有完成后验认证的结果标记为 strict certificate；其他结果使用恰当 scope。

## Potential Risks and Mitigations

1. **通用编译层在 `basic_function.jl` 中继续膨胀。**  
   Mitigation: 第一阶段优先最小改动；cluster 成为第二实际使用点后，只有接口已稳定且共享代码确有必要时才迁移到单一 `src/pauli_relaxation.jl`，避免预拆分。
2. **各模型 geometry 被迫塞进 chain/square 分支。**  
   Mitigation: 模型先拥有局部显式 site/bond/region builders；toric 与 Kitaev 第二次证明共同协议后才抽 unit-cell geometry。
3. **Heisenberg 默认 symmetry 污染带场/各向异性模型。**  
   Mitigation: 新 specification 不提供 Heisenberg 默认值；symmetry list 缺省为空，所有 generator 必须显式声明并通过 term-wise invariance。
4. **State-optimality support 增长导致不可编译或过 budget。**  
   Mitigation: 编译前生成完整 closure diagnostics；用户缩小显式 test basis，而不是编译器静默截断。
5. **Fourier metadata 错误导致“约化后 PSD”不等价。**  
   Mitigation: 先与未 Fourier 的小系统 compiled matrices 数值/符号比较；仅对用户声明的完整 orbit 启用；奇偶 momentum sectors 显式测试。
6. **固定 budget 比较不公平。**  
   Mitigation: 同时约束 max block、packed PSD entries、scalar moments 和 linear rows；adapted policy 不允许通过隐藏更多 blocks 获益。
7. **MG 能量闭合掩盖 observable 不确定。**  
   Mitigation: gate 同时要求 C(1)、C(2)、dimer order intervals；论文已显示 MG 能量精确但关联 interval 非零宽。`MinerU_markdown_2604.01555v1_2081746923717435392.md:664-699`
8. **二维远距离 observable 继续失配。**  
   Mitigation: adapted basis 显式加入 observable displacement/path/loop；报告 interval width，不只报告 energy。
9. **求解器许可证使结构测试不可运行。**  
   Mitigation: compile tests 完全无 optimizer；solver-backed tests 保持小规模并与 compile-only test 分离。
10. **`UInt16` site overflow。**  
    Mitigation: polynomial 构造时检查 `3*nsite ≤ typemax(UInt16)`；本路线不扩展 word storage，超限明确报错。
11. **严格认证无法覆盖新 block。**  
    Mitigation: scope 降为 solver bound；按 block 类型逐步扩展后验认证，不把旧 `qmb_data` 认证能力假定为自动通用。
12. **Triangular 模型缺少 frustration-free anchor。**  
    Mitigation: 明确把它定义为最终 finite-size ED stress stage，不虚构可解点；统一问题的“从 anchor 扰动”结论主要由前五个模型建立，triangular 用于检验迁移性。

## Alternative Approaches

1. **直接泛化 `GSB`。** 可立即复用成熟 Fourier、RDM 和 certification，但 `GSB` 当前把 Heisenberg basis、SU(2) 假设、state-optimality、compile 和 solve 深度耦合，风险最高。`src/bound_gsp.jl:1-2` `src/bound_gsp.jl:38-221` `src/bound_gsp.jl:226-732` 本计划不采用，保留其为 optimized reproduction backend。
2. **从 Ising 原型演化最小通用编译器（推荐）。** `_ising_moment_data` 已有清晰的 basis-product → canonical moment map 边界。`src/ising.jl:35-56` 代价是第一阶段没有全部 optimized Fourier 性能，但最适合建立正确、确定、可检查的共同路径。
3. **立即引入 PauliStrings.jl。** 可减少部分 Pauli 表示工作，但会引入 adapter 语义、依赖和 canonicalization 差异，妨碍首个 milestone；仅作为 Phase 8 可选输入层。
4. **预先设计通用 lattice/graph IR。** 可能减少后续重复，但在 toric、honeycomb、orthogonal-dimer、triangular 的真实共同接口尚未出现前容易过度设计；本计划采用模型内 geometry，第二使用点再抽取。
5. **自动 symmetry/basis discovery。** 可能提升易用性，却直接违反本路线的可审计显式设计，也难区分严格约化和启发式选择；不采用。

## Stage Gates

- [ ] **Gate A — Compiler Gate:** Phase 0–1 全部完成；compile/solve 分离、determinism、closure、diagnostics、Ising/Heisenberg 回归和 Table 2 结构测试通过。
- [ ] **Gate B — One-dimensional Anchor Gate:** MG 与解耦 dimer 两 anchor、period-2 symmetry、fixed-budget basis 对照通过。
- [ ] **Gate C — Generic Multi-body Gate:** cluster 三体+场、string observable、symmetry-breaking negative tests 通过；此时才决定共享文件抽取。
- [ ] **Gate D — Multi-sublattice Stabilizer Gate:** toric edge geometry、四体 stabilizers、短 loops、sector scope 通过。
- [ ] **Gate E — Bond-dependent Solvable Gate:** Kitaev 联合空间-轴 symmetry、flux、零场 benchmark 与小场扰动通过。
- [ ] **Gate F — Patterned Dimer Gate:** Shastry–Sutherland 扩大 unit cell、exact dimer product、非零扰动 fixed-budget 对照通过。
- [ ] **Gate G — Frustrated 2D Stress Gate:** triangular geometry/point group、ED 对照、至少两个耦合点和数值降级记录通过。

任何 gate 失败时，状态保持 `In Progress`，记录首个失败 specification、diagnostics 和最小复现；不得跳过并把后续 stage 标为 Completed。

## Final Definition of Done

- [ ] 显式 Pauli polynomial/Hamiltonian 容器、basis specification、relaxation specification、compiled artifact 和 diagnostics 已稳定并有公共或内部清晰边界。
- [ ] 程序只机械生成指定 algebra/products/constraints；不存在 auto symmetry 或 auto basis guessing。
- [ ] compile 与 solve 完全分离；compile-only tests 不依赖 Mosek。
- [ ] 通用 linear state-optimality、论文式 PSD state-optimality 和显式 RDM region 均通过 deterministic closure/Hermiticity tests。
- [ ] Ising endpoints 与 Heisenberg optimized `GSB` 回归通过，旧后端未被首阶段重写。
- [ ] Table 2 的 `N=100,d=4,r=1` sparse size 12,001 和最大 PSD block size 31 被复现；文档/测试明确 31 不是约束数。
- [ ] 六个模型按指定顺序通过各自 gate，且每个模型的 Hamiltonian、basis、symmetry、linear tests、PSD basis、RDM regions、observables、benchmark 与 stop condition 均由显式 specification 重放。
- [ ] 完整跨模型矩阵包含 uniform-local/operator-adapted、fixed budget、strengthening ablation、anchor/perturbation、局域/非局域 observable 和 result scope。
- [ ] 数值病态 diagnostics 完整；失败/降级配置不被静默删除。
- [ ] heuristic basis、equivalent reduction、valid strengthening 和 strict certificate 在类型、diagnostics 与结果表中始终分开。
- [ ] 未引入 RG、serialization、通用 qudit/fermion IR；PauliStrings.jl 若存在也仅为可选 adapter。
- [ ] 新共享源码文件均能指出至少两个真实使用点；否则共享逻辑保留在现有文件或模型专属边界。

## 清空上下文后的首条执行指令

> 从 Gate A 开始，只实施 Phase 0 与 Phase 1，不进入任何新模型：先读取 `src/QMBCertify.jl:1-43`、`src/basic_function.jl:1-509`、`src/ising.jl:1-132`、`src/bound_gsp.jl:1-890`、`test/runtests.jl:1-34` 和论文 `MinerU_markdown_2604.01555v1_2081746923717435392.md:189-616`；在不重写 `GSB` 的前提下，先补 Pauli 代数/轴编码回归，再于现有文件中实现最小显式 Pauli polynomial、basis/relaxation specification、任意 basis 的 `B†B` moment compiler、通用 linear/PSD state-optimality、显式小 RDM region、deterministic diagnostics 与 compile/solve 分离；保持旧 Ising API，并以 Ising endpoints、Heisenberg 小规模回归和 Table 2 的 `12,001 → max block 31` 作为 Gate A 验收。不要创建共享源码文件，除非在完成 Ising 与第二个实际调用点后能证明该抽取必要。