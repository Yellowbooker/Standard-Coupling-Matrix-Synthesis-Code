# Standard-Coupling-Matrix-Synthesis-Code
Standard Coupling Matrix Synthesis Code  
标准化耦合矩阵综合代码：包含广义切比雪夫滤波器耦合矩阵综合与矩阵旋转（折叠型、箭型、级联CTCQ型）代码以及基于广义切比雪夫滤波器的有耗滤波器耦合矩阵综合代码  
2025-11-24，更新六阶拓展盒型、盒型、M13=0/M24=0/M13=M24三型CQ等拓扑结构的耦合矩阵变换方法，添加三种子函数：Box6dExtract、BoxExtract、CQExtract  
请运行主函数：  
1. GChebyshevFilterMain.m（无耗广义切比雪夫滤波器耦合矩阵综合）  
2. LossyFilterMain.m（基于广义切比雪夫滤波器的有耗滤波器耦合矩阵综合）
3. PredistortionFilterMain.m（预失真滤波器耦合矩阵综合）  
子函数说明：  
1. General_Chebyshev.m（广义切比雪夫滤波函数的综合）
2. Cheby2EPF.m（广义切比雪夫滤波器对应特征多项式的综合）
3. EPF2Y.m（通过特征多项式E、P、F求解Y参数，无耗情况）
4. EPF2Ylossy.m（通过特征多项式E、P、F求解Y参数，有耗情况）
5. CheckUnitary.m（确认S参数是否满足无耗条件，即幺正性）
6. Y2CMtrans.m（Y参数到横向耦合矩阵）
7. Rotate.m（对耦合矩阵进行一次旋转，消去特定位置的矩阵元素）
8. to_foldedCM.m（旋转至折叠型耦合矩阵，顺次耦合为正，模值小于0.0001的元素置零）
9. CM2arrow.m（旋转至箭型耦合矩阵）
10. TriExtract.m（从箭型耦合矩阵中提取出CT结构，两个CT结构可转换为一个CQ结构）
11. CMFC_Response.m（计算耦合矩阵对应的S参数）
12. LossyCM2PCP.m（有耗滤波器横向耦合矩阵至并联双耦合谐振器耦合矩阵）
13. LossDistribution.m（有耗滤波器并联双耦合谐振器耦合矩阵的损耗重新分配，双耦合谐振器的Q值一致）
14. residue2.m
15. Predistortion.m（标准预失真滤波函数综合，所有极点实部同时平移相同距离）
16. AdaptivePredistortion.m（自适应预失真滤波函数综合，所有极点实部同时平移不同距离，距离由优化目标决定）
17. Box6dExtract.m（从箭型耦合矩阵中提取出六阶拓展盒型拓扑结构，具体方法见：关于六阶拓展盒型拓扑结构的说明.pdf）
18. BoxExtract.m（从箭型耦合矩阵中提取出盒型拓扑结构）
19. CQExtract.m（从箭型耦合矩阵中提取出三型CQ拓扑结构，‘13’型（M24=0），‘24’型（M13=0），‘equal’型（M13=M24））
