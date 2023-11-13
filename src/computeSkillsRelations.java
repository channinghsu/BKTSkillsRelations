import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.StreamTokenizer;
import java.util.HashMap;
import java.util.Map;

public class computeSkillsRelations {

    public String[] students_ = new String[2760000];// Number of instances
    public String[] skill_ = new String[2760000];
    public double[] right_ = new double[2760000];
    public int[] skillends_ = new int[1000];// Number of Skills
    public double[] lnsigma_ = new double[2760000];
    public int skillnum = -1;
    public int[] sourceSkillNum = new int[1000];// the number of current source skill
    public boolean lnminus1_estimation = false;
    public boolean bounded = true;
    public boolean L0Tbounded = false;
    public boolean cv = true;
    public int fold = 0;
    public int curstudent_ = 0;
    public boolean pstc = true;
    public Map<String, Double> top = new HashMap<>();
    public Map<String, Double> top_ptsc = new HashMap<>();

    public final double stepSize = 0.05;
    public final double minVal = 0.000001;
    public final Integer totalSteps = 1000000;

    class BKTParams {
        public double L0, G, S, T;

        public BKTParams(Double init) {
            if (init < 0) {
                this.L0 = Math.random() * top.get("L0");
                this.G = Math.random() * top.get("G");
                this.S = Math.random() * top.get("S");
                this.T = Math.random() * top.get("T");
            } else {
                this.L0 = init;
                this.G = init;
                this.S = init;
                this.T = init;
            }
        }

        public BKTParams(BKTParams copy, Boolean randStep) {
            this.L0 = copy.L0;
            this.G = copy.G;
            this.S = copy.S;
            this.T = copy.T;

            if (randStep) {
                Double randomchange = Math.random();
                Double thisStep = 2. * (Math.random() - 0.5) * stepSize;

                // Randomly change one of the BKT parameters.
                if (randomchange <= 0.25) {
                    this.L0 = Math.max(Math.min(this.L0 + thisStep, top.get("L0")), minVal);
                } else if (randomchange <= 0.5) {
                    this.T = Math.max(Math.min(this.T + thisStep, top.get("T")), minVal);
                } else if (randomchange <= 0.75) {
                    this.G = Math.max(Math.min(this.G + thisStep, top.get("G")), minVal);
                } else {
                    this.S = Math.max(Math.min(this.S + thisStep, top.get("S")), minVal);
                }
            }
        }

        public BKTParams(BKTParams copy) {
            this(copy, false);
        }
    }

    public double findGOOF(int start, int end, BKTParams params, boolean give_preds, boolean storelnsigma,
            int sourceskillnum) {
        double SSR = 0.0;
        String prevstudent = "FWORPLEJOHN";
        double prevL = 0.0;
        double likelihoodcorrect;
        double prevLgivenresult;
        double newL;
        curstudent_ = 0;
        Integer count = 0;

        for (int i = start; i <= end; i++) {
            if (!students_[i].equals(prevstudent)) {
                prevL = params.L0;
                prevstudent = students_[i];
                curstudent_++;
            }
            if ((!cv) || ((curstudent_ % 4 != fold) && !give_preds)) {
                if (lnminus1_estimation) {
                    likelihoodcorrect = prevL;
                } else {
                    likelihoodcorrect = (prevL * (1.0 - params.S)) + ((1.0 - prevL) * params.G);
                }
                SSR += (right_[i] - likelihoodcorrect) * (right_[i] - likelihoodcorrect);
                count++;
                prevLgivenresult = right_[i]
                        * ((prevL * (1.0 - params.S)) / ((prevL * (1 - params.S)) + ((1.0 - prevL) * (params.G))));
                prevLgivenresult += (1 - right_[i])
                        * ((prevL * params.S) / ((prevL * params.S) + ((1.0 - prevL) * (1.0 - params.G))));

                newL = prevLgivenresult + (1.0 - prevLgivenresult) * params.T;
                prevL = newL;
            }
            if ((!cv) || ((curstudent_ % 4 == fold) && give_preds)) {
                if (lnminus1_estimation) {
                    likelihoodcorrect = prevL;
                } else {
                    likelihoodcorrect = (prevL * (1.0 - params.S)) + ((1.0 - prevL) * params.G);
                }
                SSR += (right_[i] - likelihoodcorrect) * (right_[i] - likelihoodcorrect);
                count++;
                if (storelnsigma) {
                    lnsigma_[i] = likelihoodcorrect;
                }
//                System.out.print("Non-PSTC");
//                System.out.print(",");
//                System.out.print(fold);
//                System.out.print(",");
//                System.out.print("blank"); // In PSTC output, this column reflects the source skill
//                System.out.print(",");
//                System.out.print(skill_[i]);
//                System.out.print(",");
//                System.out.print(students_[i]);
//                System.out.print(",");
//                System.out.print(right_[i]);
//                System.out.print(",");
//                System.out.print("blank"); // In PSTC output, this columnn reflects the P(Ln) used for PSTC component
//                System.out.print(",");
//                System.out.println(likelihoodcorrect);

                prevLgivenresult = right_[i]
                        * ((prevL * (1.0 - params.S)) / ((prevL * (1 - params.S)) + ((1.0 - prevL) * (params.G))));
                prevLgivenresult += (1 - right_[i])
                        * ((prevL * params.S) / ((prevL * params.S) + ((1.0 - prevL) * (1.0 - params.G))));
                newL = prevLgivenresult + (1.0 - prevLgivenresult) * params.T;
                prevL = newL;
            }
        }
        return Math.sqrt(SSR / count);
    }

    public double findGOOFpstc(int start, int end, BKTParams_ptsc params, boolean give_preds, int sourceskill) {
        double SSR = 0.0;
        String prevstudent = "FWORPLEJOHN";
        double prevL = 0.0;
        double likelihoodcorrect;
        double prevLgivenresult;
        double newL;
        double plnstar;
        int j = 0;
        curstudent_ = 0;
        boolean newstudentflag;
        Integer count = 0;
        for (int i = start, sum = 0; i <= end; i++, sum++) {
            // System.out.println(SSR);
            // System.out.println(students_[i]);
            newstudentflag = false;

            if (!students_[i].equals(prevstudent)) {
                // System.out.println(students_[i]);
                prevL = params.L0;
                prevstudent = students_[i];
                curstudent_++;
                newstudentflag = true;
            }

            if ((!cv) || ((curstudent_ % 4 != fold) && !give_preds)) {
                if (lnminus1_estimation)
                    likelihoodcorrect = prevL;
                else {
                    likelihoodcorrect = (prevL * (1.0 - params.S)) + ((1.0 - prevL) * params.G);
                    likelihoodcorrect = prevL + ((1.0 - prevL) * params.T);
                    // must be adjusted for the number of observations (N - 1)
                    j = i - start + skillends_[sourceskill] - sourceSkillNum[sourceskill];
                    // 若当前知识点的答题序列大于源知识点的序列
                    if (sourceSkillNum[sourceskill] < sum) {
                        break;// 不再估计当前知识点的参数
                    } else if (newstudentflag) {
                        plnstar = likelihoodcorrect;
                    } else {
                        // 第二个观测起，源知识点开始影响当前知识点
                        plnstar = likelihoodcorrect + ((1.0 - likelihoodcorrect) * lnsigma_[j] * params.K);
                    }

                    likelihoodcorrect = plnstar;
                }
                SSR += (right_[i] - likelihoodcorrect) * (right_[i] - likelihoodcorrect);
                count++;

                prevLgivenresult = right_[i]
                        * ((prevL * (1.0 - params.S)) / ((prevL * (1 - params.S)) + ((1.0 - prevL) * (params.G))));
                prevLgivenresult += (1 - right_[i])
                        * ((prevL * params.S) / ((prevL * params.S) + ((1.0 - prevL) * (1.0 - params.G))));

                newL = prevLgivenresult + (1.0 - prevLgivenresult) * params.T;
                prevL = newL;
            }
            if ((!cv) || ((curstudent_ % 4 == fold) && give_preds)) {
                if (lnminus1_estimation)
                    likelihoodcorrect = prevL;
                else {
                    likelihoodcorrect = (prevL * (1.0 - params.S)) + ((1.0 - prevL) * params.G);
                    likelihoodcorrect = prevL + ((1.0 - prevL) * params.T);

                    j = i - start + skillends_[sourceskill] - sourceSkillNum[sourceskill];

                    // 如果源知识点的观测序列长度小于当前知识点的序列长度
                    if (sourceSkillNum[sourceskill] < sum) {
                        break; // 不再估计当前知识点的参数
                    } else if (newstudentflag) {
                        plnstar = likelihoodcorrect;
                    } else {
                        // 从第二次观测开始，源知识点开始影响当前知识点
                        plnstar = likelihoodcorrect + ((1.0 - likelihoodcorrect) * lnsigma_[j] * params.K);
                    }

                    likelihoodcorrect = plnstar;
                }
                SSR += (right_[i] - likelihoodcorrect) * (right_[i] - likelihoodcorrect);
                count++;

//                System.out.print("PSTC");
//                System.out.print(",");
//                System.out.print(fold);
//                System.out.print(",");
//                System.out.print(sourceskill);
//                System.out.print(",");
//                System.out.print(skill_[i]);
//                System.out.print(",");
//                System.out.print(students_[i]);
//                System.out.print(",");
//                System.out.print(right_[i]);
//                System.out.print(",");
//                if (j < 0) {
//                    System.out.print("skip");
//                } else {
//                    System.out.print(lnsigma_[j]);
//                }
//                System.out.print(",");
//                System.out.println(likelihoodcorrect);

                prevLgivenresult = right_[i]
                        * ((prevL * (1.0 - params.S)) / ((prevL * (1 - params.S)) + ((1.0 - prevL) * (params.G))));
                prevLgivenresult += (1 - right_[i])
                        * ((prevL * params.S) / ((prevL * params.S) + ((1.0 - prevL) * (1.0 - params.G))));

                newL = prevLgivenresult + (1.0 - prevLgivenresult) * params.T;
                prevL = newL;
            }
        }
        return Math.sqrt(SSR / count);
    }

    class BKTParams_ptsc {
        public double L0, G, S, T, K;

        public BKTParams_ptsc(Double init) {
            if (init < 0) {
                this.L0 = Math.random() * top_ptsc.get("L0");
                this.G = Math.random() * top_ptsc.get("G");
                this.S = Math.random() * top_ptsc.get("S");
                this.T = Math.random() * top_ptsc.get("T");
                this.K = Math.random() * top_ptsc.get("K");
            } else {
                this.L0 = init;
                this.G = init;
                this.S = init;
                this.T = init;
                this.K = init;
            }
        }

        public BKTParams_ptsc(BKTParams_ptsc copy, Boolean randStep) {
            this.L0 = copy.L0;
            this.G = copy.G;
            this.S = copy.S;
            this.T = copy.T;
            this.K = copy.K;

            if (randStep) {
                Double randomchange = Math.random();
                Double thisStep = 2. * (Math.random() - 0.5) * stepSize;

                // Randomly change one of the BKT parameters.
                if (randomchange <= 0.2) {
                    this.L0 = Math.max(Math.min(this.L0 + thisStep, top_ptsc.get("L0")), minVal);
                } else if (randomchange <= 0.4) {
                    this.T = Math.max(Math.min(this.T + thisStep, top_ptsc.get("T")), minVal);
                } else if (randomchange <= 0.6) {
                    this.G = Math.max(Math.min(this.G + thisStep, top_ptsc.get("G")), minVal);
                } else if (randomchange <= 0.8) {
                    this.S = Math.max(Math.min(this.S + thisStep, top_ptsc.get("S")), minVal);
                } else {
                    this.K = Math.max(Math.min(this.K + thisStep, top_ptsc.get("K")), minVal);
                }
            }
        }

        public BKTParams_ptsc(BKTParams_ptsc copy) {
            this(copy, false);
        }
    }

    public StreamTokenizer create_tokenizer(String infile) {

        try {
            StreamTokenizer st = new StreamTokenizer(new FileReader(infile));
            st.wordChars(95, 95);
            return st;
        } catch (FileNotFoundException fnfe) {
            fnfe.printStackTrace();
        }
        return null;
    }

    public void read_in_data(StreamTokenizer st_) {
        int actnum = 0;
        try {
            int tt = 724;
            skillnum = -1;
            String prevskill = "FLURG";

            tt = st_.nextToken();
            tt = st_.nextToken();

            tt = st_.nextToken();
            tt = st_.nextToken();
            tt = st_.nextToken();
            tt = st_.nextToken();
            tt = st_.nextToken();

            while (tt != StreamTokenizer.TT_EOF) {
                tt = st_.nextToken(); // num

                if (tt == StreamTokenizer.TT_EOF) {
                    prevskill = skill_[actnum - 1];
                    if (skillnum > -1)
                        skillends_[skillnum] = actnum - 1;
                    break;
                }

                tt = st_.nextToken(); // lesson

                tt = st_.nextToken();
                students_[actnum] = st_.sval;

                tt = st_.nextToken();
                skill_[actnum] = st_.sval;

                tt = st_.nextToken(); // cell

                tt = st_.nextToken();
                right_[actnum] = st_.nval;

                tt = st_.nextToken(); // eol

                // System.out.println(slip_[actnum]);

                actnum++;
                if (!skill_[actnum - 1].equals(prevskill)) {
                    prevskill = skill_[actnum - 1];
                    if (skillnum > -1)
                        skillends_[skillnum] = actnum - 2;
                    skillnum++;
                }

            }
        } catch (Exception e) {
            System.out.println(actnum);
            e.printStackTrace();
        }
    }

    public void fit_skill_model(int curskill, boolean howaboutpstc) {
        if (!howaboutpstc) {
            if (L0Tbounded) {
                top.put("L0", 0.85);
                top.put("T", 0.3);
            } else {
                top.put("L0", 0.999999);
                top.put("T", 0.999999);
            }

            if (bounded) {
                top.put("G", 0.3);
                top.put("S", 0.3);
            } else {
                top.put("G", 0.999999);
                top.put("S", 0.999999);
            }
            // oldParams is randomized.
            BKTParams oldParams = new BKTParams(-1.);
            BKTParams bestParams = new BKTParams(0.01);

            double oldRMSE = 1.;
            double newRMSE = 1.;

            double bestRMSE = 9999999.0;
            double prevBestRMSE = 9999999.0;

            double temp = 0.005;

            int startact = 0;
            if (curskill > 0) {
                startact = skillends_[curskill - 1] + 1;
            }
            int endact = skillends_[curskill];

            // Get the initial RMSE.
            oldRMSE = findGOOF(startact, endact, oldParams, false, true, curskill);

            for (Integer i = 0; i < totalSteps; i++) {
                // Take a random step.
                BKTParams newParams = new BKTParams(oldParams, true);

                newRMSE = findGOOF(startact, endact, oldParams, false, true, curskill);

                if (Math.random() <= Math.exp((oldRMSE - newRMSE) / temp)) { // Accept (otherwise move is rejected)
                    oldParams = new BKTParams(newParams);
                    oldRMSE = newRMSE;
                }

                if (newRMSE < bestRMSE) { // This method allows the RMSE to increase, but we're interested
                    bestParams = new BKTParams(newParams); // in the global minimum, so save the minimum values as the
                                                           // "best."
                    bestRMSE = newRMSE;
                }

                if (i % 10000 == 0 && i > 0) { // Every 10,000 steps, decrease the "temperature."
                    if (bestRMSE == prevBestRMSE)
                        break; // If the best estimate didn't change, we're done.

                    prevBestRMSE = bestRMSE;
                    temp = temp / 2.0;
                }

            }
            findGOOF(startact, endact, bestParams, true, true, curskill);

            System.out.print(skill_[startact]);
            System.out.print("\t");
            System.out.print(fold);
            System.out.print("\t");
            System.out.print("blank"); // In PSTC output this column reflects the source skill
            System.out.print("\t");
            System.out.print(curskill);
            System.out.print("\t");
            System.out.print(bestParams.L0);
            System.out.print("\t");
            System.out.print(bestParams.G);
            System.out.print("\t");
            System.out.print(bestParams.S);
            System.out.print("\t");
            System.out.print(bestParams.T);
            System.out.print("\t");
            System.out.print("blank"); // In PSTC output this column reflects the K-value
            System.out.println("\teol");
        }
        if (howaboutpstc) {
            top_ptsc.put("K", 0.999999);
            if (L0Tbounded) {
                top_ptsc.put("L0", 0.85);
                top_ptsc.put("T", 0.3);
            } else {
                top_ptsc.put("L0", 0.999999);
                top_ptsc.put("T", 0.999999);
            }

            if (bounded) {
                top_ptsc.put("G", 0.3);
                top_ptsc.put("S", 0.3);
                // top_ptsc.put("K", 0.6);
            } else {
                top_ptsc.put("G", 0.999999);
                top_ptsc.put("S", 0.999999);
            }
            for (int sourceskill = 0; sourceskill <= skillnum; sourceskill++) {
                if (sourceskill != curskill) {
                    // oldParams is randomized.
                    BKTParams_ptsc oldParams = new BKTParams_ptsc(-1.);
                    BKTParams_ptsc bestParams = new BKTParams_ptsc(0.01);

                    double oldRMSE = 1.;
                    double newRMSE = 1.;

                    double bestRMSE = 9999999.0;
                    double prevBestRMSE = 9999999.0;

                    double temp = 0.005;

                    int startact = 0;
                    if (curskill > 0) {
                        startact = skillends_[curskill - 1] + 1;
                    }
                    int endact = skillends_[curskill];

                    // Get the initial RMSE.
                    oldRMSE = findGOOFpstc(startact, endact, oldParams, false, sourceskill);

                    for (Integer i = 0; i < totalSteps; i++) {
                        // Take a random step.
                        BKTParams_ptsc newParams = new BKTParams_ptsc(oldParams, true);

                        newRMSE = findGOOFpstc(startact, endact, newParams, false, sourceskill);

                        if (Math.random() <= Math.exp((oldRMSE - newRMSE) / temp)) { // Accept (otherwise move is
                                                                                     // rejected)
                            oldParams = new BKTParams_ptsc(newParams);
                            oldRMSE = newRMSE;
                        }
                        // This method allows the RMSE to increase, but we're interested
                        if (newRMSE < bestRMSE) {
                            // in the global minimum, so save the minimum values as the "best."
                            bestParams = new BKTParams_ptsc(newParams);
                            bestRMSE = newRMSE;
                        }

                        if (i % 10000 == 0 && i > 0) { // Every 10,000 steps, decrease the "temperature."
                            if (bestRMSE == prevBestRMSE)
                                break; // If the best estimate didn't change, we're done.

                            prevBestRMSE = bestRMSE;
                            temp = temp / 2.0;
                        }
                    }

                    System.out.print("Source skill = ");
                    System.out.print(sourceskill);
                    System.out.print(", ");
                    System.out.print("Destination skill = ");
                    System.out.print(curskill);
                    System.out.println("\teol");

                    findGOOFpstc(startact, endact, bestParams, true, sourceskill);

                    System.out.print(skill_[startact]);
                    System.out.print("\t");
                    System.out.print(fold);
                    System.out.print("\t");
                    System.out.print(sourceskill);
                    System.out.print("\t");
                    System.out.print(curskill);
                    System.out.print("\t");
                    System.out.print(bestParams.L0);
                    System.out.print("\t");
                    System.out.print(bestParams.G);
                    System.out.print("\t");
                    System.out.print(bestParams.S);
                    System.out.print("\t");
                    System.out.print(bestParams.T);
                    System.out.print("\t");
                    System.out.print(bestParams.K);
                    System.out.println("\teol");
                } else {
                    System.out.print("skill_[startact]");
                    System.out.print("\t");
                    System.out.print(fold);
                    System.out.print("\t");
                    System.out.print("sourceskill");
                    System.out.print("\t");
                    System.out.print("curskill");
                    System.out.print("\t");
                    System.out.print("bestParams.L0");
                    System.out.print("\t");
                    System.out.print("bestParams.G");
                    System.out.print("\t");
                    System.out.print("bestParams.S");
                    System.out.print("\t");
                    System.out.print("bestParams.T");
                    System.out.print("\t");
                    System.out.print(0.99);
                    System.out.println("\teol");
                }
            }
        }
    }

    private void computelzerot(String infile_) {
        StreamTokenizer st_ = create_tokenizer(infile_);

        read_in_data(st_);
        // Add the number of source skill
        sourceSkillNum[0] = skillends_[0] + 1;
        int j = 1;
        for (int i = 1; i <= skillnum; i++) {
            sourceSkillNum[j++] = skillends_[i] - skillends_[i - 1];
        }

        System.out.println("skill\tL0\tG\tS\tT\teol");
        for (int curskill = 0; curskill <= skillnum; curskill++) {
            for (fold = 0; fold < 4; fold++) {
                fit_skill_model(curskill, false);
            }
        }

        if (pstc) {
            for (int curskill = 0; curskill <= skillnum; curskill++) {
                for (fold = 0; fold < 4; fold++) {
                    fit_skill_model(curskill, true);
                }
            }
        }
    }

    public static void main(String[] args) {
//         String infile_ = "data/as_ptsc.tsv";
//        String infile_ = "data/ct.tsv";
         String infile_ = "data/student-problem-middle.tsv";
        computeSkillsRelations c = new computeSkillsRelations();
        c.computelzerot(infile_);
    }
}
