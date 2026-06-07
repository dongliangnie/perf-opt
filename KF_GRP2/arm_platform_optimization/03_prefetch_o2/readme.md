ARM 优化阶段 3：软件预取优化

一、阶段定位

本阶段在数据对齐优化基础上继续进行 ARM 平台优化。
本阶段主要优化 dense 矩阵主线中的 mat_mul()，通过软件预取提前把后续要访问的数据加载到 Cache 中。

对应 ARM 平台优化表第 3 行：
    优化前：86.528937 ms
    优化后：70.332875 ms
    加速：18.72%
    final_sum：33407.632487

这是 ARM dense 矩阵主线中除 GCC -O2 外效果最明显的一步。

二、本阶段包含的 6 个文件

1. kf_linux.h
   与上一阶段接口一致。

2. kf_linux.c
   在 mat_mul() 中加入 __builtin_prefetch()，并对内层 j 循环进行 4 路展开。

3. kf_linux_test.c
   测试主函数。

4. timestamp.h
   时间戳接口声明。

5. timestamp.c
   时间戳接口实现。

6. readme.txt
   本阶段说明。

三、本阶段做了什么优化

上一阶段的 mat_mul() 已经使用 i-l-j 循环顺序：

    for i
        for l
            for j
                result_row[j] += a_value * b_row[j]

本阶段在内层 j 循环中加入软件预取：

    __builtin_prefetch(b_row + j + 16, 0, 1)
    __builtin_prefetch(result_row + j + 16, 1, 1)

同时将内层 j 循环做 4 路展开：

    result_row[j]     += a_value * b_row[j]
    result_row[j + 1] += a_value * b_row[j + 1]
    result_row[j + 2] += a_value * b_row[j + 2]
    result_row[j + 3] += a_value * b_row[j + 3]

四、软件预取是什么意思

矩阵乘法中，程序会连续访问大量数组数据。
如果 CPU 需要某个数据时才从内存加载，可能会等待较长时间。

软件预取的思想是：
    在真正使用某段数据之前，提前告诉 CPU：后面马上要用到这段数据。

这样 CPU 可以提前把数据加载进 Cache。
等程序真正计算到这部分数据时，数据可能已经在 Cache 中，等待时间就会减少。

五、__builtin_prefetch 参数解释

本阶段使用：

    __builtin_prefetch(b_row + j + 16, 0, 1)

含义：
    b_row + j + 16：提前预取后面即将访问的数据。
    0：表示预取用于读取。
    1：表示局部性提示，告诉编译器这段数据短期内会用到。

本阶段还使用：

    __builtin_prefetch(result_row + j + 16, 1, 1)

其中第二个参数 1 表示这个地址后面会被写入。

六、为什么本阶段加速明显

当前程序的主要瓶颈是矩阵乘法。
矩阵乘法不仅计算量大，内存访问也非常频繁。

前面 i-l-j 已经让访问更连续；
数据对齐让数组起始地址更规整；
软件预取进一步减少 CPU 等待内存的时间。

这几步叠加后，horizon=298 时由 86.528937 ms 降到 70.332875 ms。

七、需要注意的严谨说法

本阶段不应只说“纯软件预取”。
更准确的说法是：

    软件预取 + mat_mul 内层循环展开

因为代码中除了 __builtin_prefetch()，还对内层 j 循环做了 4 路展开。
所以本阶段的收益来自预取和内层展开的共同作用。

八、本阶段没有改变什么

本阶段没有改变：
    KF 数学公式
    dense 完整矩阵语义
    测试口径
    final_sum
    数据对齐设置
    编译选项 -O2

九、优化类型

本阶段属于：
    ARM 平台优化
    Cache 优化
    软件预取优化
    存储访问优化
    矩阵乘法热点优化

十、编译方法

    gcc -Wall -Wextra -O2 kf_linux.c kf_linux_test.c timestamp.c -lm -o kf_arm_prefetch_o2

十一、运行方法

    ./kf_arm_prefetch_o2 | tee arm_prefetch_o2_result.txt

十二、正确性检查

horizon=298 时，参考输出：
    time_ms 约为 70.332875
    final_sum 应为 33407.632487

如果 final_sum 与上一阶段一致，说明软件预取和内层展开没有改变 KF 计算结果。

十三、PPT 汇报建议

可以这样讲：

    在 mat_mul 的 i-l-j 连续访问基础上，我们进一步使用 __builtin_prefetch 预取 b_row 和 result_row 后续元素。
    这样可以减少 CPU 等待内存数据的时间。
    这个优化体现的是 ARM 平台存储系统和 Cache 层次上的程序优化。
