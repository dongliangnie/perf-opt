阶段 3：循环展开优化

一、阶段定位

本阶段在第 2 阶段 mat_mul 循环顺序优化的基础上继续进行算子代码优化。
本阶段仍然使用 -O0 编译，只对部分简单循环做 4 路循环展开。

对应算子优化表第 3 行：
    优化前：371.959041 ms
    优化后：370.698124 ms
    加速：0.34%
    final_sum：33407.632487

本阶段是算子优化中最后一个成功优化阶段。
后续 Vector Rotation 在本轮 clean 实验中没有加速，因此未放入本代码包。

二、本阶段包含的 6 个文件

1. kf_linux.h
   保留第 2 阶段接口和 KfWorkspace。

2. kf_linux.c
   保留 mat_mul 的 i-l-j 顺序，在此基础上对简单循环做 4 路展开。

3. kf_linux_test.c
   测试主函数。

4. timestamp.h
   时间戳接口声明。

5. timestamp.c
   时间戳接口实现。

6. readme.txt
   本阶段说明。

三、本阶段做了什么优化

本阶段主要对以下循环做 4 路展开：

1. mat_add()
   原来每次循环处理 1 个 double。
   优化后每次循环处理 4 个 double。

2. mat_sub()
   原来每次循环处理 1 个 double。
   优化后每次循环处理 4 个 double。

3. z_minus_x 计算循环
   原公式：
       z_minus_x[i] = z[i] - x[i]
   优化后每次循环处理 i、i+1、i+2、i+3。

4. x 更新循环
   原公式：
       x[i] += k_mult_z_minus_x[i]
   优化后每次循环更新 4 个元素。

5. final_sum 累加循环
   原来逐个累加 x[i]。
   优化后每次循环累加 4 个元素。

四、循环展开的原理

普通循环每一轮都需要做：

    1. 执行当前元素计算。
    2. i++。
    3. 判断 i 是否到达循环结束条件。
    4. 跳转回循环开头。

如果每次只处理 1 个元素，循环控制开销会出现很多次。

循环展开后，每轮处理 4 个元素，循环次数大约变为原来的 1/4。
这样可以减少：

    循环变量递增次数
    条件判断次数
    跳转次数

同时，展开后的代码更像连续的直线代码，也方便 CPU 流水线执行。

五、为什么本阶段加速不大

本阶段加速约 0.34%，不大但有效。
原因是：

    当前程序最大的瓶颈仍然是 mat_mul() 的三层矩阵乘法。
    mat_add()、mat_sub()、z_minus_x、x 更新、final_sum 这些循环虽然也会执行，但计算量远小于 mat_mul()。

所以循环展开能带来小幅提升，但不会像 mat_mul 循环顺序优化那样明显。

六、修改了哪些代码

主要修改 kf_linux.c：

1. mat_add()
   增加 4 路展开：
       result[i]、result[i+1]、result[i+2]、result[i+3]

2. mat_sub()
   增加 4 路展开。

3. kf_linux_iopointer() 中 z_minus_x 循环
   改成每次处理 4 个元素。

4. kf_linux_iopointer() 中 x 更新循环
   改成每次处理 4 个元素。

5. final_sum 循环
   改成每次累加 4 个元素。

七、本阶段没有改变什么

本阶段没有改变：
    KF 数学公式
    dense 完整矩阵语义
    mat_mul 的 i-l-j 循环顺序
    workspace 结构
    编译选项 -O0
    ARM 数据对齐
    软件预取
    NEON SIMD
    OpenMP

八、优化类型

本阶段属于：
    KF 算子代码优化
    循环优化
    循环展开优化
    减少循环控制开销

九、编译方法

    gcc -Wall -Wextra -O0 kf_linux.c kf_linux_test.c timestamp.c -lm -o kf_step3_loop_unroll_o0

十、运行方法

    ./kf_step3_loop_unroll_o0 | tee step3_loop_unroll_o0_result.txt

十一、正确性检查

horizon=298 时，参考输出：
    time_ms 约为 370.698124
    final_sum 应为 33407.632487

如果 final_sum 与第 2 阶段一致，说明循环展开没有改变 KF 计算结果。

十二、作为后续 ARM 优化基础的原因

本阶段是算子优化中最后一个成功阶段。
后续 ARM 平台优化从本阶段继续，而不是从 Vector Rotation 阶段继续，因为 Vector Rotation 在 clean 实验中没有取得加速。
