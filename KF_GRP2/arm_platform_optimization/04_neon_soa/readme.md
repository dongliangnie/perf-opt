ARM 优化阶段 4：SoA + NEON SIMD 专项实验

一、阶段定位

本阶段是 ARM 平台优化中的 SIMD 专项实验。
它不是前面 dense 矩阵主线的连续优化，而是单独设计的 ARM NEON 对比实验。

本阶段测试口径为：
    horizon = dim = 298
    steps = 10000
    一维逐元素 KF 更新

对应 ARM 平台优化表最后一行：
    优化前：AoS scalar，17.361479 ms
    优化后：SoA + NEON，12.291062 ms
    加速：29.20%
    final_sum：3024063.944582

二、本阶段包含的 6 个文件

1. kf_linux.h
   定义 KfMode、KfAosState、KfInput、KfOutput 和老师要求的三个接口函数。

2. kf_linux.c
   同时实现 AoS scalar 和 SoA + NEON 两种计算路径。

3. kf_linux_test.c
   测试主函数。每个 horizon 输出两行：AoS scalar 和 SoA + NEON。

4. timestamp.h
   时间戳接口声明。

5. timestamp.c
   时间戳接口实现。

6. readme.txt
   本阶段说明。

三、为什么本阶段要单独做

前面的 dense 矩阵主线使用的是完整矩阵运算：
    P、Q、R、K 都是 dim × dim 矩阵。

而 NEON 更适合处理连续数组上的逐元素运算。
所以本阶段单独构造了一维逐元素 KF 更新场景，用来验证 ARM NEON SIMD 对 KF 逐元素更新的加速效果。

注意：
    不能把本阶段写成“从 70.332875 ms 继续优化到 12.291062 ms”。
    因为前者是 dense 矩阵、steps=1；本阶段是一维逐元素更新、steps=10000。

正确写法是：
    SoA + NEON 是 ARM SIMD 专项实验，单独说明测试口径和对比对象。

四、本阶段对比对象

1. AoS scalar

AoS = Array of Structures。
每个状态元素使用一个结构体保存：

    typedef struct KfAosState_s {
        double x;
        double p;
        double q;
        double r;
        double k;
    } KfAosState;

内存布局可以理解为：

    state[0] = x0, p0, q0, r0, k0
    state[1] = x1, p1, q1, r1, k1
    state[2] = x2, p2, q2, r2, k2

这种布局对标量代码比较直观，但不适合 SIMD 连续加载 x0、x1 这类同类数据。

2. SoA + NEON

SoA = Structure of Arrays。
将同类数据分开连续存放：

    x: x0, x1, x2, x3, ...
    p: p0, p1, p2, p3, ...
    q: q0, q1, q2, q3, ...
    r: r0, r1, r2, r3, ...
    k: k0, k1, k2, k3, ...

这样 x[i] 和 x[i+1] 连续，p[i] 和 p[i+1] 连续。
这就适合 NEON 一次加载两个 double。

五、本阶段逐元素 KF 公式

每个维度 i 独立执行：

    p[i] = p[i] + q[i]
    k[i] = p[i] / (p[i] + r[i])
    x[i] = x[i] + k[i] * (z[i] - x[i])
    p[i] = (1.0 - k[i]) * p[i]

变量含义：

    x[i]：第 i 个状态维度的估计值。
    z[i]：第 i 个状态维度的测量值。
    p[i]：第 i 个维度上的估计误差。
    q[i]：过程噪声。
    r[i]：测量噪声。
    k[i]：第 i 个维度上的卡尔曼增益。

六、NEON 是怎么加速的

本阶段数据类型仍然使用 double。
在 AArch64 上，NEON 支持：

    float64x2_t

它一次可以存放两个 double。

标量代码每次处理一个 i：

    p[i]
    q[i]
    r[i]
    x[i]
    z[i]

NEON 代码每次处理两个 i：

    p[i], p[i+1]
    q[i], q[i+1]
    r[i], r[i+1]
    x[i], x[i+1]
    z[i], z[i+1]

核心 NEON 操作包括：

    vld1q_f64()    加载两个 double
    vaddq_f64()    两个 double 同时加法
    vsubq_f64()    两个 double 同时减法
    vmulq_f64()    两个 double 同时乘法
    vdivq_f64()    两个 double 同时除法
    vst1q_f64()    写回两个 double

七、为什么 SoA 和 NEON 要放在一起

SoA 单独不一定加速。
它的主要作用是让同类数据连续存放，为 NEON 连续加载创造条件。

所以本阶段在最终成功优化表中合并写成：

    SoA + NEON SIMD 优化

这样更准确。

八、代码中主要修改点

1. kf_linux.h
   新增 KfMode：
       KF_MODE_AOS_SCALAR
       KF_MODE_SOA_NEON

   新增 KfAosState：
       保存 AoS 标量版本的 x、p、q、r、k。

   KfInput 中同时包含：
       KfAosState *aos_states;
       double *x;
       double *p;
       double *q;
       double *r;
       double *k;

2. kf_linux.c
   新增两条核心路径：
       kf_core_aos_scalar()
       kf_core_soa_neon()

   kf_linux_ioself_profiling() 中分别运行：
       AoS scalar
       SoA + NEON

3. kf_linux_test.c
   输出格式变为：
       horizon,steps,time_ns,time_ms,platform,operator,method,status,final_sum

九、本阶段没有改变什么

本阶段没有改变：
    double 数据类型
    KF 逐元素更新公式
    随机种子 srand(12345)
    正确性检查方式 final_sum

十、优化类型

本阶段属于：
    ARM 平台优化
    数据布局优化
    NEON SIMD 向量化优化
    指令集层优化
    系统结构层次优化

十一、编译方法

推荐在 ARM Cortex-A72 / ARM Linux 上使用：

    gcc -Wall -Wextra -O2 -mcpu=cortex-a72 kf_linux.c kf_linux_test.c timestamp.c -lm -o kf_arm_neon_soa_o2

如果 -mcpu=cortex-a72 不支持，可以去掉该选项：

    gcc -Wall -Wextra -O2 kf_linux.c kf_linux_test.c timestamp.c -lm -o kf_arm_neon_soa_o2

十二、运行方法

    ./kf_arm_neon_soa_o2 | tee arm_neon_soa_o2_result.txt

十三、正确性检查

horizon=298，steps=10000 时，参考输出：

    AoS scalar：17.361479 ms，final_sum = 3024063.944582
    SoA + NEON：12.291062 ms，final_sum = 3024063.944582

同一个 horizon 下，两种方法的 final_sum 应保持一致。
如果 final_sum 一致，说明 SoA + NEON 只改变了数据布局和计算方式，没有改变 KF 计算结果。

十四、PPT 汇报建议

可以这样讲：

    在 ARM SIMD 专项实验中，我们对比了 AoS scalar 和 SoA + NEON。
    AoS 把每个状态的 x、p、q、r、k 放在一个结构体里，适合普通标量写法；
    SoA 把所有 x、p、q、r、k 分别连续存放，更适合 NEON 连续加载。
    在 AArch64 上，float64x2_t 一次可以处理两个 double，因此 KF 的逐元素更新可以向量化。
    实验结果显示，在 horizon=298、steps=10000 时，SoA + NEON 从 17.361479 ms 降到 12.291062 ms，加速约 29.20%。
