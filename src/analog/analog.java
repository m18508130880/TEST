package analog;

public class analog {
	// 20160825---- 多子系统雨水排水管网
	// 20150808new雨水管网过程模拟计算程序-芝加哥过程线
	// －－SQ LIU, TONGJI UNIVERSITY, 25 AUG 2016

	public static void main(String[] args) {
		// 管网基础数据：子系统数，最大管段数，最大节点数，最大管道起点数，路径最大管段数，最大计算次数，中间结果写文件指针
		int Nbr = 8;// 子系统数
		int NP = 116;// 最大管段数
		int NN = 117;// 最大节点数
		int Nstart = 2;// 最大管道起点数(自定)
		int Npline = 2;// 路径最大管段数(自定)
		int Nmax = 20;// 最大计算次数
		int Iprt = 0;// 中间结果写文件指针
		// 模拟时段数，芝加哥峰点时段位置（r=0.375），管道路径数，路径最大节点数，终点节点号
		int NT = 60;// 模拟时段数
		int NR = 23;// 芝加哥峰点时段位置（r=0.375）
		int Nroute = 2;// 管道路径数 (自定)
		int Nr_node = 2;// 路径最大节点数(自定)
		int Nend = 8;// 终点节点号

		int ip, it, k1, kp, in1, in2, in3;
		int i, j, k, kk, ik, jk, k00, kkk1 = 1;
		double dtnt, taa, tbb, AA, XX1, XX2;

		// XX[NT],qit[NT],sumqj[Nbr] [NT][NN],sumAj[Nbr] [NT][NN],Tnode[Nbr]
		// [NN][NN],sumTnode[Nbr] [NN][NN];
		double[] XX = new double[NT];
		double[] qit = new double[NT];
		double[][][] sumqj = new double[Nbr][NT][NN];
		double[][][] sumAj = new double[Nbr][NT][NN];
		double[][][] Tnode = new double[Nbr][NN][NN];
		double[][][] sumTnode = new double[Nbr][NN][NN];

		double hdc_min, hdc_max, fun_hd0, hda, hdb, hdc = 0, hde = 0, hdf, hdd;
		// hdcc[2],fun_hd[2],fun_hd0,hdcc0[NT][NP],hda,hdb,hdc,hde,hdf,hdd;
		double[] hdcc = new double[2];
		double[] fun_hd = new double[2];
		double[][] hdcc0 = new double[NT][NP];

		// 管段流速（m/s）
		// double vp[Nbr] [NP],slop[Nbr] [NP],Ad0,AD,RD,vp0=1.0;
		double[][] vp = new double[Nbr][NP];
		double[][] slop = new double[Nbr][NP];
		double Ad0, AD, RD, vp0 = 1.0;

		// double qpt[Nbr] [NT][NP],vpt[Nbr] [NT][NP],rid[Nbr]
		// [NT][NP],slopt[Nbr] [NT][NP],
		// Hwup[Nbr] [NT][NP],Hwdw[Nbr] [NT][NP],qjt[Nbr] [NT][NN],
		// overflow[Nbr] [NT][NN],Hw_over[Nbr]
		// [NT][NN],hdj0,ARD23,totalflow[Nbr] [NN],totalHw[Nbr] [NN];
		double[][][] qpt = new double[Nbr][NT][NP];
		double[][][] vpt = new double[Nbr][NT][NP];
		double[][][] rid = new double[Nbr][NT][NP];
		double[][][] slopt = new double[Nbr][NT][NP];
		double[][][] Hwup = new double[Nbr][NT][NP];
		double[][][] Hwdw = new double[Nbr][NT][NP];
		double[][][] qjt = new double[Nbr][NT][NN];
		double[][][] overflow = new double[Nbr][NT][NN];
		double[][][] Hw_over = new double[Nbr][NT][NN];
		double hdj0, ARD23 = 1;
		double[][] totalflow = new double[Nbr][NN];
		double[][] totalHw = new double[Nbr][NN];

		// 暴雨公式shanghai storm water formular:
		// (A1+C*lgP)/(t+b)**n---ln(N)=2.303log(N)---出口水位（m）
		double A1 = 17.53, C_storm = 0.95, tmin = 10, b_storm = 11.77, P_simu = 10, n_storm = 0.88, dt = 2.0, rc = 0.375;
		double[] Hw_end = new double[Nbr];

		// 节点汇水面积(ha)
		double Aj[][] = new double[Nbr][NN];
		// 节点汇水面积径流系数
		double[][] Acoef = new double[Nbr][NN];
		// 节点地面标高（m）[Nbr] [NN=23]
		double[][] Hj = new double[Nbr][NN];
		// 管网起始节点号和起始节点管底埋深<m>
		int[][] NJstart = new int[Nbr][Nstart];
		double[][] HJstart = new double[Nbr][Nstart];
		// 管网路径数和路径节点号(-99表示空节点)
		int[][][] Mroute = new int[Nbr][Nroute][Nr_node];
		// pipe branches-reverse order
		int[][][] Mbranch = new int[Nbr][Nstart][Npline];
		// 管段上游节点号I0,下游节点号J0，管段长度(m),摩|阻系数
		int[][] I0 = new int[Nbr][NP];
		int[][] J0 = new int[Nbr][NP];
		double[][] lp = new double[Nbr][NP];
		// double slp[Nbr] [NP]={};
		double[][] slp = new double[Nbr][NP];
		// 管段直径(m)，上游管底高程(m)，下游管底高程(m)
		double[][] dpl = new double[Nbr][NP];
		double[][] ZJup = new double[Nbr][NP];
		double[][] ZJdw = new double[Nbr][NP];
		// double[] ZJdw = new double[Nbr];
		// ----------------------------------------------------------------------------------------------------------
		//
		// --输出数据文件开始---
		/*
		 * ofstream outfile; outfile.open("20160808-雨水管网过程模拟-芝加哥过程线-10-60.txt");
		 */
		// ================= 赋初值 ===============

		// ofstream outfile;
		// outfile.open("20160808-雨水管网过程模拟-芝加哥过程线-10-60.txt");
		System.out.println("20160808-雨水管网过程模拟-芝加哥过程线-10-60.txt");
		// ================= 赋初值 ===============
		for (kk = 0; kk < Nbr; kk++) {
			for (i = 0; i < NT; i++) {
				for (j = 0; j < NN; j++)
					sumAj[kk][i][j] = 0;
			}
			for (i = 0; i < NT; i++) {
				for (j = 0; j < NN; j++)
					sumqj[kk][i][j] = 0;
			}
			for (i = 0; i < NN; i++) {
				for (j = 0; j < NN; j++) {
					if (i == j) {
						Tnode[kk][i][j] = 0;
					} else {
						Tnode[kk][i][j] = -99;
					}
				}
			}
			for (i = 0; i < NN; i++) {
				for (j = 0; j < NN; j++)
					sumTnode[kk][i][j] = 0;
			}
		}
		// ==================Tnode-sumTnode=========================
		for (kk = 0; kk < Nbr; kk++) {
			for (i = 0; i < NP; i++)
				vp[kk][i] = vp0;
			for (kp = 0; kp < NP; kp++) {
				in1 = I0[kk][kp];
				in2 = J0[kk][kp];
				Tnode[kk][in1][in2] = lp[kk][kp] / vp[kk][kp] / 60;
				slop[kk][kp] = (ZJup[kk][kp] - ZJdw[kk][kp]) / lp[kk][kp];
			}
			//
			for (i = 0; i < Nroute; i++) {
				for (j = 0; j < Nr_node; j++) {
					in1 = Mroute[kk][i][j];
					if (in1 >= 0) {
						for (k = j + 1; k < Nr_node; k++) {
							in2 = Mroute[kk][i][k - 1];
							in3 = Mroute[kk][i][k];
							if (in3 >= 0) {
								sumTnode[kk][in1][in3] = sumTnode[kk][in1][in2]
										+ Tnode[kk][in2][in3];
							}
						}
					}
				}
			}
		}
		// kk
		// =====print Mroute[i][j], Tnode, sumTnode,Mbranch[i][j]====
		// cout<<"pipe no.  I0    J0"<<endl;
		System.out.println("pipe no.  I0    J0");
		for (kk = 0; kk < Nbr; kk++) {
			for (i = 0; i < NP; i++) {
				System.out.println("i,I0[kk][i],J0[kk][i]:" + i + I0[kk][i]
						+ J0[kk][i]);
				// s.Format("%6d%6d%6d\n", i,I0[kk] [i],J0[kk] [i]);
				// cout<<s;
			}
			// outfile<<endl;
			//
			// outfile<<" ip=";
			System.out.println("ip=");
			for (i = 0; i < NP; i++) {
				System.out.println(i);
				// s.Format("%4d",i);
				// outfile<<s;
			}
			// outfile<<endl;
			//
			// outfile<<" I0=";
			System.out.println(" I0=");
			for (i = 0; i < NP; i++) {
				System.out.println("I0[kk] [i]" + I0[kk][i]);
				// s.Format("%4d",I0[kk] [i]);

				// outfile<<s;
			}
			// outfile<<endl;
			//
			// outfile<<" J0=";
			System.out.println(" J0=");
			for (i = 0; i < NP; i++) {
				System.out.println("J0[kk] [i]" + J0[kk][i]);
				// s.Format("%4d",J0[kk] [i]);
				// outfile<<s;
			}
			// outfile<<endl;
		}// kk
			//
			// outfile<<endl;
			// outfile<<"===========  print Mroute[kk] [i][j]"<<endl;
		System.out.println("===========  print ");
		for (kk = 0; kk < Nbr; kk++) {
			for (i = 0; i < Nroute; i++) {
				for (j = 0; j < Nr_node; j++) {
					System.out.println("Mroute[kk][i][j]" + Mroute[kk][i][j]);
					// s.Format("%6d",Mroute[kk] [i][j]);
					// outfile<<s;
				}
				// outfile<<endl; }
			}// kk
				//
				// outfile<<endl;
				// outfile<<"===========  print Mbranch[kk] [i][j]"<<endl;
			System.out.println("===========  print Mbranch[kk] [i][j]");
			for (kk = 0; i < Nbr; kk++) {
				for (i = 0; i < Nstart; i++) {
					for (j = 0; j < Npline; j++) {
						System.out.println("Mbranch[kk] [i][j]"
								+ Mbranch[kk][i][j]);
						// s.Format("%6d",Mbranch[kk] [i][j]);
					}
					// outfile<<endl; }
				}// kk
					//
					// outfile<<"===========  print Tnode[kk] [i][j]"<<endl;
					// outfile<<"====j=  "<<endl;
					// outfile<<"      ";
				System.out.println("===========  print Tnode[kk] [i][j]");
				System.out.println("====j=  ");
				for (kk = 0; i < Nbr; kk++) {
					for (j = 0; j < NN; j++) { // s.Format("%6d",j);
												// outfile<<s;
						System.out.println(j);
					}
					// outfile<<endl;
					//
					for (i = 0; i < NN; i++) {
						if (i < 10) {
							// outfile<<"i="<<i<<"   ";
							System.out.println("i=" + i);
						} else {
							// outfile<<"i="<<i<<"  ";
							System.out.println("i=" + i);

						}
						for (j = 0; j < NN; j++) {
							if (Tnode[kk][i][j] < 0.0) {
								// outfile<<"      ";
								// System.out.println("        ");
							} else {
								// s.Format("%6.2f",Tnode[kk] [i][j]);
								System.out.println("Tnode[kk] [i][j]"
										+ Tnode[kk][i][j]);
								// outfile<<s;
							}
						}
						// outfile<<endl;
					}
				}// kk
					//
					// outfile<<endl;
					// outfile<<"===========  print sumTnode[kk] [i][j]"<<endl;
					// outfile<<"==j=  ";
					// outfile<<"      ";
				System.out.println("===========  print sumTnode[kk] [i][j]");
				System.out.println("==j=  ");
				for (kk = 0; i < Nbr; kk++) {
					for (j = 0; j < NN; j++) { // s.Format("%6d",j);
												// outfile<<s;
						System.out.println("j=" + j);
					}
					// outfile<<endl;
					//
					for (i = 0; i < NN; i++) {// outfile<<"i="<<i<<"   ";

						System.out.println("i=" + i);
						for (j = 0; j < NN; j++) {
							if (sumTnode[kk][i][j] <= 0.0) {// outfile<<"      ";

								// System.out.println("      ");
							} else

							{ // s.Format("%6.2f",sumTnode[kk] [i][j]);
								// outfile<<s;
								System.out.println("sumTnode[kk] [i][j]"
										+ sumTnode[kk][i][j]);

							}
						}
						// outfile<<endl;
					}
				}// kk
					// ================= 管网稳态流动模拟===================
					//
					// -----------动态模拟流量计算------------------------------------
					// ----------------节点汇水面积(ha)和汇水流量(m3/sec)计算--------
					// outfile<<endl;
					// outfile<<"===========  管网动态模拟计算      重现期＝ "<<P_simu<<"  年   时段数＝ "<<NT<<"       终点水位＝ "<<Hw_end<<"  m  ========="<<endl;
				System.out.println("===========  管网动态模拟计算      重现期＝ " + P_simu
						+ "  年   时段数＝ " + NT + "       终点水位＝ " + Hw_end
						+ "  m  =========");
				// ----------rainfall intensity at every time step-----------
				// xxxxxxx
				// 芝加哥过程线
				AA = A1 + A1 * C_storm * Math.log(P_simu) / 2.303;
				for (it = 0; it < NT; it++) {
					if (it <= NR) {
						dtnt = dt * (float) (it);
						tbb = dt * (float) (NR) - dtnt;
						XX1 = AA * ((1.0 - n_storm) * tbb / rc + b_storm);
						XX2 = Math.pow((tbb / rc + b_storm), (n_storm + 1.0));
					} else {
						dtnt = dt * (float) (it);
						taa = dtnt - dt * (float) (NR);
						XX1 = AA
								* ((1.0 - n_storm) * taa / (1.0 - rc) + b_storm);
						XX2 = Math.pow((taa / (1.0 - rc) + b_storm),
								(n_storm + 1.0));
					}
					XX[it] = XX1 / XX2;
					qit[it] = 167.0 * XX[it] / 1000.0;
				}
				//
				// outfile<<endl;
				// outfile<<"    it      dtnt      XX[it]     qit[it]"<<endl;
				System.out.println("\t\tit    \tdtnt\t\tXX[it]\t\tqit[it]");
				for (it = 0; it < NT; it++) {
					dtnt = dt * (float) (it);
					// s.Format("%6d%10.2lf%12.6lf%12.6lf\n",it,dtnt,XX[it],qit[it]);
					// outfile<<s;
					System.out.println("\t\t" + it + "\t    " + dtnt + "\t\t"
							+ XX[it] + "\t" + qit[it]);
				}
				// outfile<<endl;
				// xxxxxxx
				//
				for (kk = 0; i < Nbr; kk++) {
					for (it = 0; it < NT; it++) {
						dtnt = dt + dt * (float) (it);
						for (j = 0; j < NN; j++) {
							sumAj[kk][it][j] = Aj[kk][j];
							sumqj[kk][it][j] = Aj[kk][j] * qit[it]
									* Acoef[kk][j];
							for (i = 0; i < NN; i++) {
								if (sumTnode[kk][i][j] > 0
										&& sumTnode[kk][i][j] < dtnt) {
									sumAj[kk][it][j] = sumAj[kk][it][j]
											+ Aj[kk][i];
									sumqj[kk][it][j] = sumqj[kk][it][j]
											+ Aj[kk][i] * qit[it]
											* Acoef[kk][i];
								}
							}
						}
					}
					// print sumAj[it][j] and sumqj[it][j]
					// outfile<<"  sumAj[kk] [it][j]="<<endl;
					System.out.println("    sumAj[kk] [it][j]=");
					for (it = 0; it < NT; it++) {
						for (j = 0; j < NN; j++) {// s.Format("%8.2lf",sumAj[kk]
													// [it][j]);
													// outfile<<s;
							System.out.println(sumAj[kk][it][j]);
						}
						// outfile<<endl;
					}
					// outfile<<endl;
					//
					// outfile<<"  sumqj[kk] [it][j]="<<endl;
					System.out.println("    sumqj[kk] [it][j]=");
					for (it = 0; it < NT; it++) {
						for (j = 0; j < NN; j++) {// s.Format("%8.2lf",sumqj[kk]
													// [it][j]);
													// outfile<<s;
							System.out.println(sumqj[kk][it][j]);
						}
						// outfile<<endl;
					}
					// outfile<<endl;
				}// kk
					// ---------------------------------------------------------------
				for (it = 0; it < NT; it++) {
					for (i = 0; i < NN; i++) {
						overflow[kk][it][i] = 0.0;
						Hw_over[kk][it][i] = 0.0;
					}
				}
				for (it = 0; it < NT; it++) {
					for (j = 0; j < NP; j++)
						qpt[kk][it][j] = -99.0;
				}
				// ---------------------------------------------------------------
				// kk-----000开始
				for (kk = 0; i < Nbr; kk++) {
					for (it = 0; it < NT; it++)
					// --1
					{// outfile<<" it="<<it<<"  qpt[kk] [it][k]=";
						System.out.println("   it=" + it
								+ "    qpt[kk] [it][k]=");
						for (j = 0; j < NN; j++) {
							for (k = 0; k < NP; k++) {
								if (I0[kk][k] == j) {
									qpt[kk][it][k] = sumqj[kk][it][j];
									// s.Format("%8.2lf",qpt[kk] [it][k]);
									// outfile<<s;
									System.out.println(qpt[kk][it][k]);
								}
							}
						}
						// outfile<<endl;
						// -------------------???????????????090127?????????????????------------------------
						for (ik = 0; ik < Nstart; ik++)
						// --2
						{
							for (jk = 0; jk < Npline; jk++)
							// --3
							{
								kp = Mbranch[kk][ik][jk];
								if (kp >= 0)
								// --4
								{
									if (J0[kk][kp] == Nend) {
										Hwdw[kk][it][kp] = Hw_end[kk];
										// outfile<<"   it= "<<it<<"   kp= "<<kp<<"  Hwdm= "<<Hwdw[kk][it][kp]<<"  Hw_end= "<<Hw_end<<endl;
										System.out.println("   it= " + it
												+ "   kp= " + kp + "  Hwdm= "
												+ Hwdw[kk][it][kp]
												+ "  Hw_end= " + Hw_end);
									} else {
										for (k1 = 0; k1 < NP; k1++) {
											if (I0[kk][k1] == J0[kk][kp])
												Hwdw[kk][it][kp] = Hwup[kk][it][k1];
										}
									}
									//
									Ad0 = 0.7854 * Math.pow(dpl[kk][kp], 2.0);
									hdj0 = ZJdw[kk][kp] + dpl[kk][kp];
									if (Hwdw[kk][it][kp] >= hdj0) {
										hdcc0[it][kp] = 1.0;
										rid[kk][it][kp] = dpl[kk][kp] / 4.0;
										vpt[kk][it][kp] = qpt[kk][it][kp] / Ad0;
										slopt[kk][it][kp] = 10.29
												* Math.pow(slp[kk][kp], 2.0)
												* Math.pow(qpt[kk][it][kp], 2.0)
												/ Math.pow(dpl[kk][kp], 5.333);
										Hwup[kk][it][kp] = Hwdw[kk][it][kp]
												+ slopt[kk][it][kp]
												* lp[kk][kp];
										if (Hwup[kk][it][kp] >= Hj[kk][kp]) {
											Hwup[kk][it][kp] = Hj[kk][kp];
											slopt[kk][it][kp] = (Hwup[kk][it][kp] - Hwdw[kk][it][kp])
													/ lp[kk][kp];
											vpt[kk][it][kp] = Math.pow(
													rid[it][kk][kp], 0.6667)
													* Math.pow(
															slopt[kk][it][kp],
															0.5) / slp[kk][kp];
											qpt[kk][it][kp] = vpt[kk][it][kp]
													* Ad0;
										}
									} else
									// --5
									{
										hdc_min = (Hwdw[kk][it][kp] - ZJdw[kk][kp])
												/ dpl[kk][kp];
										slopt[kk][it][kp] = slop[kk][kp];
										if (hdc_min < 0.0)
											hdc_min = 0.0;
										ARD23 = slp[kk][kp]
												* qpt[kk][it][kp]
												/ Math.pow(slopt[kk][it][kp],
														0.5);
										// ----------h/d---0.618aaa---------------------------
										hdc_max = 1.0;
										fun_hd0 = 1000.0;
										for (kkk1 = 1; kkk1 < Nmax
												&& fun_hd0 > 0.0001; kkk1++)
										// --6
										{
											if (kkk1 == 1) {
												hdcc[0] = hdc_min;
												hdcc[1] = 1.0;
											} else {
												hdcc[0] = hdc_min
														+ (hdc_max - hdc_min)
														* 0.382;
												hdcc[1] = hdc_min
														+ (hdc_max - hdc_min)
														* 0.618;
											}
											for (k00 = 0; k00 < 2; k00++) {
												hda = 1.0 - 2 * hdcc[k00];
												hdb = 1.0 - hdcc[k00];
												hdc = Math.pow(Math.cos(hda),
														-1);
												hdd = hdcc[k00] * hdb;
												hde = hda * Math.pow(hdd, 0.5);
												hdf = Math
														.pow(dpl[kk][kp], 2.0);
												AD = hdf / 4.0 * hdc - hdf
														/ 2.00 * hde;
												RD = dpl[kk][kp] / 4.0
														- dpl[kk][kp] * hde
														/ (2.0 * hdc);
												fun_hd[k00] = AD
														* Math.pow(RD,
																2.0 / 3.0)
														- ARD23;
											}
											if (Math.abs(fun_hd[0]) > Math
													.abs(fun_hd[1])) {
												hdc_min = hdcc[0];
												fun_hd0 = Math.abs(fun_hd[0]);
												hdcc0[it][kp] = hdcc[1];
											} else {
												hdc_max = hdcc[1];
												fun_hd0 = Math.abs(fun_hd[1]);
												hdcc0[it][kp] = hdcc[0];
											}
											//
											if (Iprt == 1) {// outfile<<"it=  "<<it<<"  ik=  "<<ik<<"jk=  "<<jk<<"  kp=  "<<kp<<"  kkk1=  "<<kkk1<<"  hdcc[0]=  "<<hdcc[0]<<"  hdcc[1]=  "<<hdcc[1]
												// <<"  fun_hd[0]= "<<fun_hd[0]<<"  fun_hd[1]= "<<fun_hd[1]<<"  ARD23= "<<ARD23<<endl;
												System.out.println("it=  " + it
														+ "  ik=  " + ik
														+ "jk=  " + jk
														+ "  kp=  " + kp
														+ "  kkk1=  " + kkk1
														+ "  hdcc[0]=  "
														+ hdcc[0]
														+ "  hdcc[1]=  "
														+ hdcc[1]
														+ "  fun_hd[0]= "
														+ fun_hd[0]
														+ "  fun_hd[1]= "
														+ fun_hd[1]
														+ "  ARD23= " + ARD23);
											}
										}// 6-kkk1
											//
										rid[kk][it][kp] = (0.25 - 0.5 * hde
												/ hdc)
												* dpl[kk][kp];
										Hwdw[kk][it][kp] = ZJdw[kk][kp]
												+ hdcc0[it][kp] * dpl[kk][kp];
										vpt[kk][it][kp] = Math.pow(
												rid[kk][it][kp], 2.0 / 3.0)
												* Math.pow(slopt[kk][it][kp],
														0.5) / slp[kk][kp];
										Hwup[kk][it][kp] = Hwdw[kk][it][kp]
												+ slopt[kk][it][kp]
												* lp[kk][kp];
										// ---------h/d---0.618aaa---------------------------
									} // 5--end_hdcc0[it][kp]
										//
										// outfile<<"   it= "<<it<<"   kp= "<<kp<<"   I0[kp]= "<<I0[kk]
										// [kp]<<"  Hwdm= "<<Hwdw[kk][it][kp]<<"  Hj= "<<Hj[I0[kk]
										// [kp]]
										// <<"  hdcc0= "<<hdcc0[it][kp]<<"  qpt= "<<qpt[kk][it][kp]<<"  vpt= "<<vpt[kk][it][kp]<<"  hdcc0= "<<hdcc0[it][kp]<<" fun_hd0= "<<fun_hd0<<endl;
									System.out.println("it=  " + it + "  ik=  "
											+ ik + "jk=  " + jk + "  kp=  "
											+ kp + "  kkk1=  " + kkk1
											+ "  hdcc[0]=  " + hdcc[0]
											+ "  hdcc[1]=  " + hdcc[1]
											+ "  fun_hd[0]= " + fun_hd[0]
											+ "  fun_hd[1]= " + fun_hd[1]
											+ "  ARD23= " + ARD23);
								}// 4 if(kp>=0) end
							}// 3 ---jk end ---
						}// --2---ik end ---Hj[kp]
							// --------------------------------- 开始计算溢流节点
							// ---------------
							//
						for (i = 0; i < NP; i++) {
							j = I0[kk][i];
							if (Hwup[kk][it][i] == Hj[kk][j]) {
								overflow[kk][it][j] = overflow[kk][it][j]
										- qpt[kk][it][i];
								for (ip = 0; ip < NP; ip++) {
									k1 = J0[kk][ip];
									if (k1 == j)
										overflow[kk][it][j] = overflow[kk][it][j]
												+ qpt[kk][it][ip];
								}
								qjt[kk][it][j] = Aj[kk][j] * qit[it]
										* Acoef[kk][j];
								overflow[kk][it][j] = (overflow[kk][it][j] + qjt[kk][it][j])
										* dt * 60.0;
								Hw_over[kk][it][j] = overflow[kk][it][j]
										/ Aj[kk][j] / 10000.0 * 1000.0;
							}
						}
						// ----qjt[it][j],overflow[it][j],totalflow[j],totalHw[j]
						// ------------------ 计算溢流节点结束 ---------------
						// outfile<<endl;
						// outfile<<"    it   管段号  I0   J0   管径dpl   管段qp 水力半径R    充满度¨¨ 流速¨′(m/s)  上游水位  下游水位  上管底高  下管底高  上地面高"<<endl;
						System.out
								.println("    it   管段号  I0   J0   管径dpl   管段qp 水力半径R    充满度¨¨ 流速¨′(m/s)  上游水位  下游水位  上管底高  下管底高  上地面高");
						for (i = 0; i < NP; i++) {// s.Format("%6d%6d%6d%5d%10.2lf%10.3lf%10.3lf%10.3lf%10.3lf%10.3lf%10.3lf%10.3lf%10.3lf%10.3lf\n",
													// it,i,I0[kk] [i],J0[kk]
													// [i],dpl[kk] [i],qpt[kk]
													// [it][i],rid[kk]
													// [it][i],hdcc0
													// [it][i],vpt[kk]
													// [it][i],Hwup[kk]
													// [it][i],Hwdw[kk]
													// [it][i],ZJup[kk]
													// [i],ZJdw[kk] [i],Hj[kk]
													// [I0[i]]);
							// outfile<<s;
							System.out.println(it + i + I0[kk][i] + J0[kk][i]
									+ dpl[kk][i] + qpt[kk][it][i]
									+ rid[kk][it][i] + hdcc0[it][i]
									+ vpt[kk][it][i] + Hwup[kk][it][i]
									+ Hwdw[kk][it][i] + ZJup[kk][i]
									+ ZJdw[kk][i] + Hj[kk][I0[kk][i]]);
						}
						// outfile<<endl;
					}// 1-- it end ---
						// ------- 计算节点积水量和积水深度¨¨(m)----
					for (j = 0; j < NN; j++) {
						totalflow[kk][j] = 0.0;
						totalHw[kk][j] = 0.0;
					}
					for (j = 0; j < NN; j++) {
						for (it = 0; it < NT; it++) {
							totalflow[kk][j] = totalflow[kk][j]
									+ overflow[kk][it][j];
							totalHw[kk][j] = totalHw[kk][j]
									+ Hw_over[kk][it][j];
						}
					}
					// -----屏幕输出管段水力计算结束------
					// cout<<" it 管段号 I0 J0 管径dpl 管段qp 水力半径R 充满度 流速(m/s) 上游水位
					// 下游水位 上管底高 下管底高 上地面高”<<endl;
					System.out
							.println("    it   管段号  I0   J0   管径dpl   管段qp 水力半径R    充满度 流速(m/s)  上游水位  下游水位  上管底高  下管底高  上地面高");
					// outfile<<endl;
					for (it = 0; it < NT; it++) {// outfile<<"  it= "<<it<<endl;
													// outfile<<endl;
													// outfile<<" it 管段号 I0 J0
													// 管径dpl 管段qp 水力半径R 充满度¨¨
													// 流速(m/s) 上游水位 下游水位 上管底高程
													// 下管底高程<<endl;
						for (i = 0; i < NP; i++) {// s.Format("%6d%6d%6d%5d%10.2lf%10.3lf%10.3lf%10.3lf%10.3lf%10.3lf%10.3lf%10.3lf%10.3lf%10.3lf\n",
													// it,i,I0[kk] [i],J0[kk]
													// [i],dpl[kk] [i],qpt[kk]
													// [it][i],rid[kk]
													// [it][i],hdcc0
													// [it][i],vpt[kk]
													// [it][i],Hwup[kk]
													// [it][i],Hwdw[kk]
													// [it][i],ZJup[kk]
													// [i],ZJdw[kk] [i],Hj[kk]
													// [I0[i]]);
													// outfile<<s;
							// cout<<s;
							System.out.println(it + i + I0[kk][i] + J0[kk][i]
									+ dpl[kk][i] + qpt[kk][it][i]
									+ rid[kk][it][i] + hdcc0[it][i]
									+ vpt[kk][it][i] + Hwup[kk][it][i]
									+ Hwdw[kk][it][i] + ZJup[kk][i]
									+ ZJdw[kk][i] + Hj[kk][I0[kk][i]]);
						}
						// outfile<<endl;
					}
					// outfile<<endl;
					// ------------- 节点溢流计算结果输出 ---------------
					// outfile<<" ======== 时段节点积水量(m3) ========"<<endl;
					// outfile<<"  i=    ";
					System.out
							.println(" ======== 时段节点积水量(m3) ========\n i=    ");
					for (i = 0; i < NN; i++) {
						if (i < 10) {
							// outfile<<"  "<<i<<"   ";
							System.out.print("  " + i + "  ");
						} else {
							// outfile<<" "<<i<<"   ";
							System.out.print("  " + i + "  ");
						}
					}
					// outfile<<endl;
					// outfile<<"it="<<endl;
					System.out.println("it=");
					for (it = 0; it < NT; it++) {
						if (it < 10) {
							// outfile<<" "<<it<<"   ";
							System.out.print("  " + it + "  ");
						} else {
							// outfile<<it<<"   ";
							System.out.print("  ");
						}
						//
						for (i = 0; i < NN; i++) {
							if (overflow[kk][it][i] <= 0.0) {
								// outfile<<"      ";
								System.out.print("  ");
							} else {
								// s.Format("%6.1f",overflow[kk] [it][i]);
								// outfile<<s;
								System.out.println(overflow[kk][it][i]);
							}
						}
						// outfile<<endl;
					}
					//
					// outfile<<" ======== 时段节点积水深度输出(mm) ========"<<endl;
					// outfile<<"  i=    ";
					System.out.print(" ======== 时段节点积水深度输出(mm) ========\ni=  ");
					for (i = 0; i < NN; i++) {
						if (i < 10) {
							// outfile<<"  "<<i<<"   ";
							System.out.print("  " + i + "  ");
						} else {
							// outfile<<" "<<i<<"   ";
							System.out.print("  " + i + "  ");
						}
					} // outfile<<endl;
						// outfile<<"it="<<endl;
					System.out.println("it=");
					for (it = 0; it < NT; it++) {
						if (it < 10) {
							// outfile<<" "<<it<<"   ";
							System.out.print(" " + it + "   ");
						} else {
							// outfile<<it<<"   ";
							System.out.print(" " + it + "   ");
						}
						//
						for (i = 0; i < NN; i++) {
							if (overflow[kk][it][i] <= 0.0) {
								// outfile<<"      ";
								System.out.print("  ");
							} else {
								// s.Format("%6.1f",Hw_over[kk] [it][i]);
								// outfile<<s;
								System.out.println(Hw_over[kk][it][i]);
							}
						}
						// outfile<<endl;
					}
				}// kk---000结束
					//
					// outfile.close();

			}
		}
	}
}
