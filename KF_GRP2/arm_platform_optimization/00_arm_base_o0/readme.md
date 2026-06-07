ARM 优化阶段 0：ARM 平台优化基准版本

一、阶段定位

本阶段是 ARM 平台优化部分的起点。
代码来自算子优化阶段的最后一个成功版本：
    operator_optimization/03_loop_unroll

也就是说，本阶段已经包含算子层面的成功优化：
    workspace 预分配
    mat_mul 循环顺序优化
    循环展开

但是本阶段不新增任何 ARM 平台优化。
本阶段仍然使用 -O0 编译，作为 ARM 平台优化表的 baseline。

参考结果：
    horizon=298
    time_ms：370.698124 ms
    final_sum：33407.632487

二、本阶段包含的 6 个文件

1. kf_linux.h
   KF 对外接口和 KfWorkspace 声明。

2. kf_linux.c
   来自算子优化最终成功版本，包含 workspace、mat_mul i-l-j、循环展开。

3. kf_linux_test.c
   测试主函数。

4. timestamp.h
   时间戳接口声明。

5. timestamp.c
   时间戳接口实现。

6. readme.txt
   本阶段说明。

三、为什么 ARM 优化从这个阶段开始

我们把实验分为两大部分：

第一部分是算子优化：
    先不启用编译器优化，使用 -O0 看 C 代码本身优化带来的效果。

第二部分是 ARM 平台优化：
    在算子代码已经优化到比较合理的基础上，再使用 ARM 平台相关优化。

因此，ARM 平台优化不从老师基础版直接开始，而是从算子优化最终成功版开始。
这样逻辑更清楚：
    先优化算子代码，再利用 ARM 平台特性继续优化。

四、本阶段已有的算子优化

1. workspace 预分配
   工作矩阵由 KfWorkspace 统一管理，不在核心函数中反复 calloc/free。

2. mat_mul 循环顺序优化
   mat_mul 使用 i-l-j 顺序，让 b_row 和 result_row 连续访问。

3. 循环展开
   mat_add、mat_sub、z_minus_x、x 更新、final_sum 等循环做 4 路展开。

五、本阶段没有做哪些 ARM 优化

本阶段没有做：
    GCC -O2
    数据对齐
    软件预取
    NEON SIMD
    OpenMP
    分支预测

六、优化类型

本阶段不算 ARM 优化，只是 ARM 优化的基准版本。

七、编译方法

    gcc -Wall -Wextra -O0 kf_linux.c kf_linux_test.c timestamp.c -lm -o kf_arm_base_o0

八、运行方法

    ./kf_arm_base_o0 | tee arm_base_o0_result.txt

九、正确性检查

horizon=298 时，参考输出：
    time_ms 约为 370.698124
    final_sum 应为 33407.632487

后续 ARM dense 主线阶段都应保持 final_sum = 33407.632487。
