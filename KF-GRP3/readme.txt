卡尔曼滤波器多版本优化实现
一、项目概述
-----------
本项目针对线性卡尔曼滤波器（对角协方差假设）实现了多个优化版本，
用于对比不同算法、数据类型、向量化指令和并行策略对计算性能的影响。

所有版本均保持相同的外部接口（kf_linux.h），通过编译链接不同的
实现文件（.c）即可得到不同性能特性的可执行程序。

二、文件说明
-----------
1. 接口与通用模块
   kf_linux.h              外部接口声明（KfInput, KfOutput, RetCode）
   timestamp.h / .c        高精度计时工具（基于 clock_gettime）
   kf_linux_test.c         通用测试程序，循环调用性能测试接口

2. 原始矩阵版本（基线）
   kf_linux.c              全矩阵运算，O(n^3) 复杂度，double 精度
                           （包含矩阵加法、乘法、求逆等通用操作）

3. 对角优化版本（标量）
   kf_linux_diag.c         将对角协方差退化为向量运算，O(n) 复杂度
                           仍保留部分动态分配，double 精度
   kf_linux_diag_scalar.c  标量极致优化版：合并循环、常量折叠、
                           无中间数组，double 精度
   kf_linux_float_diag_scalar.c  同上，但数据类型为 float，进一步提升
                                  吞吐量并降低内存占用

4. NEON 向量化版本（ARM 平台）
   kf_linux_neon.c         NEON 向量化矩阵乘法（double），循环分块
   kf_linux_neon_omp.c     NEON 向量化 + OpenMP 多核并行（double）
   kf_linux_float_neon_omp.c   float + NEON + OpenMP 全并行版本
                              （注意：此版本为全矩阵运算，非对角优化）

三、编译方法
-----------
系统要求：
  - Linux 操作系统
  - GCC 编译器（支持 C99、OpenMP、ARM NEON 内部函数）
  - ARM 平台（对于 NEON 版本）或 x86 平台（基线/标量版本）
  - 若编译 NEON 版本，需确保 CPU 支持 NEON 扩展

1. 编译原始矩阵版本（x86/ARM 均可）
   gcc -O3 -march=native -o test_matrix \
       kf_linux.c kf_linux_test.c timestamp.c -lm

2. 编译对角 double 版本（标量）
   gcc -O3 -march=native -o test_diag \
       kf_linux_diag_scalar.c kf_linux_test.c timestamp.c -lm

3. 编译对角 float 版本（标量）
   gcc -O3 -march=native -o test_float_diag \
       kf_linux_float_diag_scalar.c kf_linux_test.c timestamp.c -lm
   注意：此处测试程序仍会调用 kf_linux_ioself_profiling，内部将 float
         结果写回 double 输出，精度会有所损失。

4. 编译 NEON 版本（仅 ARM，如树莓派、Jetson 等）
   # 纯 NEON（无 OpenMP）
   gcc -O3 -march=armv8-a+simd -o test_neon \
       kf_linux_neon.c kf_linux_test.c timestamp.c -lm

   # NEON + OpenMP（双精度）
   gcc -O3 -march=armv8-a+simd -fopenmp -o test_neon_omp \
       kf_linux_neon_omp.c kf_linux_test.c timestamp.c -lm

   # NEON + OpenMP（单精度 float）
   gcc -O3 -march=armv8-a+simd -fopenmp -o test_float_neon_omp \
       kf_linux_float_neon_omp.c kf_linux_test.c timestamp.c -lm

   对于 ARMv7（32位），可将 -march=armv8-a+simd 替换为 -mfpu=neon。

四、运行测试
-----------
直接运行生成的可执行文件，将输出不同维度下的单次滤波耗时（纳秒）。
示例：
   ./test_matrix
   dim=10, time=xxxx ns
   dim=19, time=xxxx ns
   ...
   （维度从 10 开始，每次增加 9，直到 300）

若使用 OpenMP 版本，可通过环境变量控制线程数：
   OMP_NUM_THREADS=4 ./test_neon_omp

五、重要注意事项
---------------
1. 对角假设：标量优化版本（diag*）严格依赖协方差矩阵为对角阵的假设，
   若输入非对角矩阵将导致错误结果。原始矩阵版本和 NEON 版本无此限制。

2. 内存对齐：NEON 向量加载指令（vld1q_f64 / vld1q_f32）要求数据地址
   16 字节对齐。若使用 malloc 分配，建议替换为 aligned_alloc 或
   posix_memalign 以保证对齐，否则可能触发总线错误或性能下降。

3. 精度与类型：
   - float 版本（kf_linux_float_diag_scalar, kf_linux_float_neon_omp）
     将内部计算全部改为 float，输出时转换为 double。长时间迭代可能
     积累舍入误差，建议对比 double 版本验证一致性。
   - 在 NEON + OpenMP 的 float 版本中，因接口仍使用 KfOutput (double*)，
     内部通过强制类型转换实现，调用者必须确保传入的实际内存为 float 数组，
     否则会导致内存越界。

4. 健壮性宏：源码中包含若干条件编译的健壮性保护（如除零保护、协方差非负
   钳位、NaN/Inf 检查）。如需启用，请在编译时添加 -DNUMERIC_ROBUST 或
   -DENABLE_NAN_CHECK 等宏定义。

5. 伪随机数：为便于比较，所有版本均使用固定种子（srand(12345)）生成
   模拟测量值。在实际应用中，测量值应由传感器提供，移除 rand() 调用。

6. 性能测试：kf_linux_test.c 中的循环仅测试单次滤波耗时，未包含内存分配
   时间。若需评测整体开销，可修改测试程序。

7. 平台兼容性：timestamp.c 使用 POSIX 时钟接口（clock_gettime），
   适用于 Linux。若移植到其他系统，需替换时间测量实现。

六、参考资料
-----------
项目代码注释中包含了详细的优化思路与健壮性建议，可作为学习和扩展的参考。