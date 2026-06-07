KF 卡尔曼滤波算子 ARM Linux 程序性能优化代码包

一、代码包总体说明

本代码包用于课程研讨题目：
    KF 卡尔曼滤波算子在 ARM Linux 平台上的程序性能优化

本代码包严格按照老师要求整理，每一个阶段都是独立文件夹，每个文件夹中都包含老师要求的 6 个文件：
    1. kf_linux.h
    2. kf_linux.c
    3. kf_linux_test.c
    4. timestamp.h
    5. timestamp.c
    6. readme.txt

每个阶段都可以单独进入文件夹后编译运行，不依赖其他文件夹中的源码。

二、为什么分成两个大文件夹

本实验最终分成两条主线：

1. operator_optimization
   基于 KF 算子本身的代码优化。
   这一部分主要研究：在不依赖 ARM 特殊硬件优化的情况下，如何通过改造 C 代码本身提高 KF 算子速度。
   统一使用 -O0 编译，目的是尽量减少编译器自动优化的影响，看清楚“代码结构本身”带来的优化效果。

2. arm_platform_optimization
   基于 ARM 平台和系统结构的优化。
   这一部分主要研究：在已经完成算子代码优化后，如何利用 ARM Linux 平台上的编译器、Cache、数据对齐、软件预取、NEON SIMD 等系统结构特性进一步提升性能。

三、operator_optimization 文件夹说明

该文件夹包含 4 个阶段：

0. 00_base_required
   老师要求的 6 文件基础版。
   只完成接口规范和文件结构，不做性能优化。
   这是算子优化的干净 baseline。

1. 01_workspace_prealloc
   workspace 预分配优化。
   将 kf_linux_iopointer() 内部反复申请的工作数组整理到 KfWorkspace 中，由外部提前分配和初始化。

2. 02_matmul_cache
   mat_mul 矩阵乘法循环顺序优化。
   将矩阵乘法循环顺序由 i-j-l 改为 i-l-j，改善 C 语言行优先存储下的连续访存效果。

3. 03_loop_unroll
   循环展开优化。
   对 mat_add、mat_sub、z_minus_x、x 更新、final_sum 等简单循环做 4 路展开。

说明：
    Vector Rotation 双缓冲区轮换我们做过实验，但本轮 clean 实验中未取得加速，因此不放入成功优化代码包。

四、operator_optimization 推荐运行顺序

建议依次运行：

    cd operator_optimization/00_base_required
    gcc -Wall -Wextra -O0 kf_linux.c kf_linux_test.c timestamp.c -lm -o kf_step0_o0
    ./kf_step0_o0 | tee step0_base_o0_result.txt

    cd ../01_workspace_prealloc
    gcc -Wall -Wextra -O0 kf_linux.c kf_linux_test.c timestamp.c -lm -o kf_step1_workspace_o0
    ./kf_step1_workspace_o0 | tee step1_workspace_o0_result.txt

    cd ../02_matmul_cache
    gcc -Wall -Wextra -O0 kf_linux.c kf_linux_test.c timestamp.c -lm -o kf_step2_matmul_cache_o0
    ./kf_step2_matmul_cache_o0 | tee step2_matmul_cache_o0_result.txt

    cd ../03_loop_unroll
    gcc -Wall -Wextra -O0 kf_linux.c kf_linux_test.c timestamp.c -lm -o kf_step3_loop_unroll_o0
    ./kf_step3_loop_unroll_o0 | tee step3_loop_unroll_o0_result.txt

五、算子优化参考性能表

测试口径：
    platform = ARM Linux
    horizon = dim = 298
    steps = 1
    编译选项 = -O0
    dense 完整矩阵 KF

参考结果如下：

    Step 0 基础版：526.086937 ms，final_sum = 33407.632487
    Step 1 workspace 预分配：524.148291 ms，final_sum = 33407.632487
    Step 2 mat_mul 循环顺序优化：371.959041 ms，final_sum = 33407.632487
    Step 3 循环展开优化：370.698124 ms，final_sum = 33407.632487

成功优化表：

    1. workspace 预分配
       526.086937 ms -> 524.148291 ms
       加速约 0.37%

    2. mat_mul 循环顺序优化
       524.148291 ms -> 371.959041 ms
       加速约 29.04%

    3. 循环展开优化
       371.959041 ms -> 370.698124 ms
       加速约 0.34%

算子优化整体效果：
    526.086937 ms -> 370.698124 ms
    整体加速约 29.54%

六、arm_platform_optimization 文件夹说明

该文件夹包含 5 个阶段：

0. 00_arm_base_o0
   ARM 平台优化基准版本。
   代码来自算子优化最终成功版本，即 03_loop_unroll。
   使用 -O0 编译，作为 ARM 优化的起点。

1. 01_compiler_o2
   GCC -O2 编译器优化。
   不改源码，只将编译选项从 -O0 改为 -O2。

2. 02_alignment_o2
   数据对齐优化。
   使用 posix_memalign() 将高频访问数组按 64 字节对齐。

3. 03_prefetch_o2
   软件预取优化。
   在 mat_mul() 中加入 __builtin_prefetch()，提前预取后续要访问的 b_row 和 result_row。

4. 04_neon_soa
   SoA + NEON SIMD 专项实验。
   本阶段不是 dense 矩阵主线的连续优化，而是 ARM SIMD 专项测试。
   对比 AoS scalar 和 SoA + NEON，一次处理两个 double。

七、arm_platform_optimization 推荐运行顺序

ARM dense 主线：

    cd arm_platform_optimization/00_arm_base_o0
    gcc -Wall -Wextra -O0 kf_linux.c kf_linux_test.c timestamp.c -lm -o kf_arm_base_o0
    ./kf_arm_base_o0 | tee arm_base_o0_result.txt

    cd ../01_compiler_o2
    gcc -Wall -Wextra -O2 kf_linux.c kf_linux_test.c timestamp.c -lm -o kf_arm_o2
    ./kf_arm_o2 | tee arm_compiler_o2_result.txt

    cd ../02_alignment_o2
    gcc -Wall -Wextra -O2 kf_linux.c kf_linux_test.c timestamp.c -lm -o kf_arm_alignment_o2
    ./kf_arm_alignment_o2 | tee arm_alignment_o2_result.txt

    cd ../03_prefetch_o2
    gcc -Wall -Wextra -O2 kf_linux.c kf_linux_test.c timestamp.c -lm -o kf_arm_prefetch_o2
    ./kf_arm_prefetch_o2 | tee arm_prefetch_o2_result.txt

ARM SIMD 专项：

    cd ../04_neon_soa
    gcc -Wall -Wextra -O2 -mcpu=cortex-a72 kf_linux.c kf_linux_test.c timestamp.c -lm -o kf_arm_neon_soa_o2
    ./kf_arm_neon_soa_o2 | tee arm_neon_soa_o2_result.txt

如果平台不支持 -mcpu=cortex-a72，可以先去掉该选项编译，但在 ARM AArch64 平台上建议保留。

八、ARM 平台优化参考性能表

A. dense 矩阵主线

测试口径：
    platform = ARM Linux
    horizon = dim = 298
    steps = 1
    dense 完整矩阵 KF

参考结果：

    00_arm_base_o0：370.698124 ms，final_sum = 33407.632487
    01_compiler_o2：89.478562 ms，final_sum = 33407.632487
    02_alignment_o2：86.528937 ms，final_sum = 33407.632487
    03_prefetch_o2：70.332875 ms，final_sum = 33407.632487

成功优化表：

    1. GCC -O2 编译器优化
       370.698124 ms -> 89.478562 ms
       加速约 75.86%

    2. 数据对齐优化
       89.478562 ms -> 86.528937 ms
       加速约 3.30%

    3. 软件预取 + mat_mul 内层展开
       86.528937 ms -> 70.332875 ms
       加速约 18.72%

ARM dense 主线整体效果：
    370.698124 ms -> 70.332875 ms
    整体加速约 81.03%

B. SoA + NEON SIMD 专项

测试口径：
    platform = ARM Linux
    horizon = dim = 298
    steps = 10000
    一维逐元素 KF 更新

参考结果：
    AoS scalar：17.361479 ms，final_sum = 3024063.944582
    SoA + NEON：12.291062 ms，final_sum = 3024063.944582
    加速约 29.20%

注意：
    SoA + NEON 阶段与 dense 矩阵主线测试口径不同，不能写成“从 70.332875 ms 继续优化到 12.291062 ms”。
    正确写法是：SoA + NEON 是 ARM SIMD 专项实验，单独说明测试口径和对比对象。

九、正确性检查方法

每个阶段输出中都有 final_sum。
同一测试口径下，如果 final_sum 保持一致，说明优化没有改变 KF 计算结果。

对于 dense 矩阵主线：
    horizon=298 时 final_sum 应为 33407.632487。

对于 SoA + NEON 专项：
    horizon=298，steps=10000 时 final_sum 应为 3024063.944582。



