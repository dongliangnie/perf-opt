ARM 优化阶段 2：数据对齐优化

一、阶段定位

本阶段在 GCC -O2 编译器优化基础上继续进行 ARM 平台优化。
本阶段主要利用 ARM 存储系统和 Cache 访问特点，对高频访问数组进行 64 字节对齐。

对应 ARM 平台优化表第 2 行：
    优化前：89.478562 ms
    优化后：86.528937 ms
    加速：3.30%
    final_sum：33407.632487

二、本阶段包含的 6 个文件

1. kf_linux.h
   与上一阶段接口一致。

2. kf_linux.c
   新增 aligned_calloc()，将主要数组改为 64 字节对齐分配。

3. kf_linux_test.c
   测试主函数。

4. timestamp.h
   时间戳接口声明。

5. timestamp.c
   时间戳接口实现。

6. readme.txt
   本阶段说明。

三、本阶段做了什么优化

上一阶段中，数组使用普通 calloc() 分配。
本阶段改为：

    aligned_calloc()
        内部调用 posix_memalign(&ptr, 64, total_size)
        再使用 memset(ptr, 0, total_size) 清零

也就是说，本阶段把高频访问的数据按 64 字节边界对齐。

涉及数组包括：

    measurements
    x
    p
    q
    r
    k
    identity
    temp_p_plus_r
    inv_temp_p_plus_r
    temp_i_minus_k
    p_old
    z_minus_x
    k_mult_z_minus_x

四、为什么要做数据对齐

ARM CPU 访问内存时，不是每次只读取一个 double，而是以 Cache line 为单位把一段连续数据加载到 Cache 中。
如果数据起始地址不规整，可能出现：

    1. 一次访问跨越更多 Cache line。
    2. 连续访问效率下降。
    3. 后续向量化或连续加载不够友好。

64 字节对齐可以让数组起始地址更规整，连续访问时更容易被 Cache 高效加载。

五、为什么本阶段属于 ARM 平台优化

数据对齐本质上是根据处理器存储系统特点进行优化。
它关注的是：

    Cache line
    内存访问对齐
    连续加载效率
    存储系统性能

这些都是系统结构层次的内容，因此本阶段归入 ARM 平台优化，而不是单纯的算法优化。

六、修改了哪些代码

主要修改 kf_linux.c：

1. 文件开头加入：

    #define _POSIX_C_SOURCE 200112L

2. 新增宏：

    #define KF_ALIGNMENT_BYTES 64

3. 新增函数：

    static void *aligned_calloc(size_t count, size_t size)

4. 将普通 calloc 改为 aligned_calloc，例如：

    input->measurements = aligned_calloc(...)
    workspace->p = aligned_calloc(...)
    workspace->q = aligned_calloc(...)
    workspace->r = aligned_calloc(...)

七、本阶段没有改变什么

本阶段没有改变：
    KF 数学公式
    mat_mul 循环顺序
    循环展开逻辑
    测试口径
    编译选项 -O2
    dense 完整矩阵计算语义

八、优化类型

本阶段属于：
    ARM 平台优化
    数据对齐优化
    Cache 访问优化
    存储系统优化

九、编译方法

    gcc -Wall -Wextra -O2 kf_linux.c kf_linux_test.c timestamp.c -lm -o kf_arm_alignment_o2

十、运行方法

    ./kf_arm_alignment_o2 | tee arm_alignment_o2_result.txt

十一、正确性检查

horizon=298 时，参考输出：
    time_ms 约为 86.528937
    final_sum 应为 33407.632487

如果 final_sum 与上一阶段一致，说明数据对齐只改变了内存分配方式，没有改变 KF 计算结果。
