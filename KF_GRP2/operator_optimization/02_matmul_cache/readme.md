阶段 2：mat_mul 矩阵乘法循环顺序优化

一、阶段定位

本阶段在第 1 阶段 workspace 预分配基础上继续进行算子代码优化。
本阶段仍然使用 -O0 编译，只改变 mat_mul() 的循环顺序，不做其他优化。

这是算子优化阶段中最关键的一步。
因为 dense 完整矩阵 KF 的主要瓶颈就在矩阵乘法。

对应算子优化表第 2 行：
    优化前：524.148291 ms
    优化后：371.959041 ms
    加速：29.04%
    final_sum：33407.632487

二、本阶段包含的 6 个文件

1. kf_linux.h
   保留第 1 阶段的接口和 KfWorkspace 声明。

2. kf_linux.c
   只修改 mat_mul()，将循环顺序由 i-j-l 改为 i-l-j。

3. kf_linux_test.c
   测试主函数。

4. timestamp.h
   时间戳接口声明。

5. timestamp.c
   时间戳接口实现。

6. readme.txt
   本阶段说明。

三、原始 mat_mul() 为什么慢

第 1 阶段中的 mat_mul() 使用的是 i-j-l 循环顺序：

    for i
        for j
            for l
                result[i][j] += a[i][l] * b[l][j]

在 C 语言中，二维矩阵实际存储为一维数组，并且是按行连续存储。
也就是说：

    b[l * m + j]

当内层循环 l 变化、j 固定时，程序访问的是 B 矩阵的同一列。
但是 C 语言中同一列的元素在内存中不是连续的。
这会导致：

    1. 访存跳跃大。
    2. Cache 命中率低。
    3. CPU 经常等待内存数据。
    4. 矩阵规模越大，影响越明显。

当 horizon=298 时，mat_mul() 的计算量非常大，因此这个问题会成为主要瓶颈。

四、本阶段做了什么优化

本阶段把 mat_mul() 改为 i-l-j 循环顺序：

    for i
        for l
            for j
                result[i][j] += a[i][l] * b[l][j]

优化后，在最内层 j 循环中：

    b_row[j] 连续访问
    result_row[j] 连续访问

这更符合 C 语言按行存储的内存布局。

五、为什么这样不改变计算结果

矩阵乘法的数学公式是：

    C[i][j] = sum(A[i][l] * B[l][j])

原来的 i-j-l 和优化后的 i-l-j，本质上只是改变了循环遍历顺序。
每个 C[i][j] 仍然累加相同的 A[i][l] * B[l][j]。
因此数学结果不变。

本阶段 final_sum 仍然是 33407.632487，也证明了结果没有改变。

六、修改了哪些代码

主要修改 kf_linux.c 中的 mat_mul()。

原来写法：

    memset(result, 0, sizeof(double) * n * m);
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < m; j++) {
            double sum_value = 0.0;
            for (size_t l = 0; l < k_common; l++) {
                sum_value += a[i * k_common + l] * b[l * m + j];
            }
            result[i * m + j] = sum_value;
        }
    }

优化后写法：

    memset(result, 0, sizeof(double) * n * m);
    for (size_t i = 0; i < n; i++) {
        double *result_row = result + i * m;
        for (size_t l = 0; l < k_common; l++) {
            double a_value = a[i * k_common + l];
            const double *b_row = b + l * m;
            for (size_t j = 0; j < m; j++) {
                result_row[j] += a_value * b_row[j];
            }
        }
    }

七、为什么本阶段加速明显

本阶段没有减少矩阵乘法的理论乘加次数，复杂度仍然是 O(N^3)。
但是它明显改善了内存访问方式。

CPU 不是只看计算次数，还要看数据能不能及时从内存和 Cache 中取到。
优化后连续访问 b_row 和 result_row，空间局部性更好，所以速度明显提升。

八、本阶段没有改变什么

本阶段没有改变：
    KF 数学公式
    dense 完整矩阵语义
    workspace 结构
    mat_add()、mat_sub() 循环方式
    循环展开
    Vector Rotation
    编译选项 -O0
    ARM 平台专项优化

九、优化类型

本阶段属于：
    KF 算子代码优化
    矩阵乘法实现优化
    Cache 访问顺序优化
    数据局部性优化

十、编译方法

    gcc -Wall -Wextra -O0 kf_linux.c kf_linux_test.c timestamp.c -lm -o kf_step2_matmul_cache_o0

十一、运行方法

    ./kf_step2_matmul_cache_o0 | tee step2_matmul_cache_o0_result.txt

十二、正确性检查

horizon=298 时，参考输出：
    time_ms 约为 371.959041
    final_sum 应为 33407.632487

如果 final_sum 与第 1 阶段一致，说明循环顺序优化没有改变 KF 计算结果。
