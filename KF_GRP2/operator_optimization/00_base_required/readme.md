阶段 0：老师要求的 6 文件基础版

一、阶段定位

本阶段是整个实验的起点，也是后续所有算子优化的干净 baseline。
本阶段只完成老师要求的代码接口和文件结构，不做任何性能优化。

也就是说，本阶段的目标不是让程序更快，而是先把老师给的单文件 KF 程序整理成符合研讨要求的工程形式。
后续所有优化都应该在这个基础上一步一步修改，避免不同优化混在一起导致实验数据不清楚。

本阶段对应算子优化表中的第 0 阶段：
    编译选项：-O0
    测试口径：horizon=dim=298，steps=1，dense 完整矩阵 KF
    参考时间：526.086937 ms
    final_sum：33407.632487

二、本阶段包含的 6 个文件

1. kf_linux.h
   对外接口声明文件。
   里面包含老师要求的 RetCode_t、KfInput、KfOutput，以及三个接口函数声明。

2. kf_linux.c
   KF 算子实现文件。
   本阶段保留 dense 完整矩阵计算逻辑，工作矩阵仍然在 kf_linux_iopointer() 内部申请和释放。

3. kf_linux_test.c
   测试主函数。
   main() 从 horizon=10 开始，每次增加 9，一直测试到 KF_MAX_DIM=300 附近，即最大 horizon=298。

4. timestamp.h
   时间戳接口声明。
   提供 TimeStamp、timestamp()、timestamp_diff()。

5. timestamp.c
   时间戳接口实现。
   使用 clock_gettime(CLOCK_MONOTONIC, ...) 获取纳秒级计时。

6. readme.txt
   当前阶段说明文件。

三、三个接口函数说明

1. RetCode_t kf_linux_iopointer(int horizon, void *input, void *output)

   这是核心 KF 算子函数。
   调用者传入 input 和 output 指针，函数内部检查参数是否合法，然后执行 KF 计算。
   本阶段中，该函数内部会申请 x、P、Q、R、K、临时矩阵等工作空间，计算完成后再释放。

2. RetCode_t kf_linux_ioself_profiling(int horizon)

   这是带计时的自测接口。
   它自己生成输入数据，然后调用 kf_linux_iopointer()，只对 kf_linux_iopointer() 的执行过程计时，最后打印结果。
   输出格式为：
       horizon,time_ns,time_ms,platform,operator,status,final_sum

3. RetCode_t kf_linux_ioself(int horizon)

   这是功能测试接口。
   它自己生成输入并调用核心函数，但不打印时间，适合只检查功能是否正确。

四、horizon 在本实验中的含义

本实验中：
    horizon = dim

也就是说，horizon 表示 KF 状态维度，也表示矩阵规模。
例如 horizon=298 时，P、Q、R、K 等矩阵都是 298 × 298。

本阶段 steps_as_iterations = 1。
也就是说，每个 horizon 只执行一次 KF 更新。

五、本阶段实现的 KF 数学模型

老师原始代码是一个简化 KF：
    F = I
    H = I
    状态维度 = 测量维度
    P、Q、R 初始化为对角阵

代码中保留的核心公式是：

    P = P + Q
    K = P * inv(P + R)
    x = x + K * (z - x)
    P = (I - K) * P

变量含义：

    x：当前状态估计值。
    z：当前测量值。
    P：估计误差协方差，表示当前估计的不确定性。
    Q：过程噪声协方差，表示系统预测过程中的不确定性。
    R：测量噪声协方差，表示传感器测量的不确定性。
    K：卡尔曼增益，决定更相信预测值还是测量值。
    I：单位矩阵。

简单理解：
    先用 P = P + Q 增加预测不确定性；
    再计算 K，让程序判断测量值和估计值各占多大权重；
    然后用 x = x + K * (z - x) 修正状态估计；
    最后用 P = (I - K) * P 更新误差协方差。

六、本阶段代码特点

1. 保留完整矩阵存储。
   P、Q、R、K、identity、临时矩阵都按 dim × dim 分配。

2. 保留普通矩阵乘法。
   mat_mul() 使用原始 i-j-l 三层循环。

3. 保留核心函数内部申请和释放工作空间。
   kf_linux_iopointer() 内部有多次 calloc/free。

4. 保留老师原始简化 KF 逻辑。
   不为了性能改写数学公式。

七、本阶段没有做哪些优化

本阶段没有做：

    workspace 预分配
    mat_mul 循环顺序优化
    循环展开
    Vector Rotation
    GCC -O2 编译器优化
    数据对齐
    软件预取
    NEON SIMD
    OpenMP

这很重要，因为它保证第 0 阶段是干净的 baseline。

八、优化类型

本阶段不算性能优化。
本阶段属于：
    工程规范化
    接口标准化
    老师要求的 6 文件结构整理

九、编译方法

本阶段作为算子优化 baseline，统一使用 -O0：

    gcc -Wall -Wextra -O0 kf_linux.c kf_linux_test.c timestamp.c -lm -o kf_step0_o0

十、运行方法

    ./kf_step0_o0 | tee step0_base_o0_result.txt

十一、正确性检查

参考结果：
    horizon=298
    time_ms 约为 526.086937
    final_sum 应为 33407.632487

不同机器上 time_ms 会波动，但同一套输入下 final_sum 应保持一致。
如果后续优化后的 final_sum 仍然是 33407.632487，说明优化没有改变 KF 计算结果。
