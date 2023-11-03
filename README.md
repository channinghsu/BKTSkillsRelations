# BKTSkillsRelations

本项目在两个现有项目的基础上进行开发，分别是 `BKTSimulatedAnnealing` 和 `BKT-PSTC`。它利用 `BKTSimulatedAnnealing` 中的模拟退火算法来估计 BKT 参数，并应用 `BKT-PSTC` 中的方法来计算知识点之间的关系。最终的输出包括每个知识点的 BKT 参数以及它们之间的关系因子 K。

## 改进之处

1. 提高算法效率：用模拟退火算法替代 `BKT-PSTC` 中的暴力搜索算法，大幅减少计算时间。

2. 适用于不同数据集：在 `BKT-PSTC` 中，原本要求在估计参数 K 时，数据集中所有知识点的答题数必须相同。这意味着如果想计算不同数据集中知识点之间的关系，必须先将所有知识点的答题数调整为相同的数量，这一过程相当繁琐。在本项目的改进中，不再对答题数的要求如此严格。如果当前知识点的答题数大于源知识点的答题数，我们将不再估计当前知识点的参数（T、G、S、L、K），从而使得输入不同答题数的知识点变得更加方便。

   原始代码：

   ```java
   j = i - start + sourceskillend - 804 - 1; // j 为源知识点的第（当前观测 - 1）个观测，其中804 - 1为源知识点的观测序列长度
   
   if (newstudentflag) 
       plnstar = likelihoodcorrect;
   else {
       // 从第二次观测开始，源知识点开始影响当前知识点
       plnstar = likelihoodcorrect + ((1.0 - likelihoodcorrect) * lnsigma_[j] * K);
   }
   ```

   改进后的代码：

   ```java
      // j 为源知识点的第（当前观测 - 1）个观测
      j = i - start + skillends_[sourceskill] - sourceSkillNum[sourceskill];
   
      if (sourceSkillNum[sourceskill] < sum) { // 如果源知识点的观测序列长度小于当前知识点的序列长度
          break; // 不再估计当前知识点的参数
      } else if (newstudentflag) {
          plnstar = likelihoodcorrect;
      } else {
          // 从第二次观测开始，源知识点开始影响当前知识点
          plnstar = likelihoodcorrect + ((1.0 - likelihoodcorrect) * lnsigma_[j] * params.K);
      }
   ```

   