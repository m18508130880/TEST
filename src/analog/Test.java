package analog;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.text.DecimalFormat;

public class Test {
	// 20160825---- 多子系统雨水排水管网
	// 20150808new雨水管网过程模拟计算程序-芝加哥过程线
	// －－SQ LIU, TONGJI UNIVERSITY, 25 AUG 2016

	public static void main(String[] args) throws FileNotFoundException {
		
		String FileName = "201608026-雨水管网过程模拟-芝加哥过程线-华家池演示002.txt";
		FileOutputStream fs = new FileOutputStream(new File(FileName));
		PrintStream printStream = new PrintStream(fs);
		printStream.println(FileName);
		
		DecimalFormat df=new DecimalFormat("##.##");

		// 管网基础数据：子系统数，最大管段数，最大节点数，最大管道起点数，路径最大管段数，最大计算次数，中间结果写文件指针
		// int Nbr = 8;//子系统数
		int NP = 9;// 管段数
		int NN = 10;// 节点数
		int Nstart = 3;// 管道起点数(自定)
		int Npline = 7;// 路径最大管段数(自定)
		int Nmax = 20;// 最大计算次数
		int Iprt = 0;// 中间结果写文件指针
		// 模拟时段数，芝加哥峰点时段位置（r=0.375），管道路径数，路径最大节点数，终点节点号
		int NT = 60;// 模拟时段数
		int NR = 23;// 芝加哥峰点时段位置（r=0.375）
		int Nroute = 3;// 管道路径数 (自定)
		int Nr_node = 8;// 路径最大节点数(自定)
		int Nend = 8;// 终点节点号
		
/** 管井数据加载	
 *  汇水面积Aj	
 *  径流系数Acoef 
 *  地面标高 Hj
 *  井底标高Hj
 				顶部标高	底部标高	直径					
	YJ002001 	5.244 	3.894 	0.7 	YG002000,  		YG002001  		起  点  	铁混  	原始探测 	华家池校区演示  	 
	YJ002002 	5.191 	3.842 	0.7 	YG002001,  		YG002002  		中间点  	铁混  	原始探测 	华家池校区演示  	 
	YJ002003 	5.177 	3.784 	0.7 	YG002002,  		YG002003  		中间点  	铁混  	原始探测 	华家池校区演示  	 
	YJ002004 	5.208 	3.733 	0.7 	YG002003,  		YG002004  		中间点  	铁混  	原始探测 	华家池校区演示  	 
	YJ002005 	5.221 	3.687 	0.7 	YG002004,  		YG002005  		中间点  	铁混  	原始探测 	华家池校区演示  	 
	YJ002006 	5.201 	3.643 	0.7 	YG002005,  		YG002006  		中间点  	铁混  	原始探测 	华家池校区演示  	 
	YJ002007 	5.2 	3.25 	0.7 	YG002006,YG002009 YG002007  	中间点  	铁混  	原始探测 	华家池校区演示  	 
	YJ002008 	5.121 	3.171 	0.7 	YG002007,YG002008  YG002999  	终  点  	铁混  	原始探测 	华家池校区演示  	 
	YJ002009 	5.131 	3.731 	0.7 	YG002000,  		YG002008  		起  点  	铁混  	原始探测 	华家池校区演示  	 
	YJ002010 	5.186 	3.486 	0.7 	YG002000,  		YG002009  		起  点  	铁混  	原始探测 	华家池校区演示  	
*/
		// 管网起始节点号和起始节点管底埋深<m>
//		int[] NJstart = new int[]{1, 9, 10};
//		double[] HJstart = new double[]{3.894, 3.731, 3.486};		
		// 节点汇水面积(ha)3.5
		double Aj[] = new double[]{3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5};
		// 节点汇水面积径流系数0.6
		double[] Acoef = new double[]{0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6};
		// 节点地面标高（m）[NN=23]
		double[] Hj = new double[]{5.244,5.191,5.177,5.208,5.221,5.201,5.2,5.121,5.131,5.186};
		
/** 管段数据加载
 *  起点管井I0
 *  终点管井J0
 *  长度lp
 *  直径dpl
 *  摩阻系数slp
 *  起端标高ZJup
 *  终端标高ZJdw
				直径		长度								起端标高	终端标高
	YG002001 	0.3  	28.5  	YJ002001 	YJ002002  	3.894 	3.842 	PE  	 	原始探测 	华家池校区演示  	 
	YG002002 	0.3  	32  	YJ002002 	YJ002003  	3.842 	3.784 	PE  	 	原始探测 	华家池校区演示  	 
	YG002003 	0.3  	28.6  	YJ002003 	YJ002004  	3.784 	3.733 	PE  	 	原始探测 	华家池校区演示  	 
	YG002004 	0.3  	25.4  	YJ002004 	YJ002005  	3.733 	3.687 	PE  	 	原始探测 	华家池校区演示  	 
	YG002005 	0.3  	24.7  	YJ002005 	YJ002006  	3.687 	3.643 	PE  	 	原始探测 	华家池校区演示  	 
	YG002006 	0.3  	23.5  	YJ002006 	YJ002007  	3.643 	3.601 	PE  	 	原始探测 	华家池校区演示  	 
	YG002007 	0.3  	30.4  	YJ002007 	YJ002008  	3.601 	3.546 	PE  	 	原始探测 	华家池校区演示  	 
	YG002008 	0.3  	15.5  	YJ002009 	YJ002008  	3.731 	3.171 	PE  	 	原始探测 	华家池校区演示  	 
	YG002009 	0.3  	4.3  	YJ002010 	YJ002007  	3.886 	3.7 	PE  	 	原始探测 	华家池校区演示  	
 */
		// 管段上游节点号I0,下游节点号J0,
		int[] I0 = new int[]{1, 2, 3, 4, 5, 6, 7, 9, 10};
		int[] J0 = new int[]{2, 3, 4, 5, 6, 7, 8, 8, 7};
		//管段长度
		double[] lp = new double[]{28.5, 32, 28.6, 25.4, 24.7, 23.5, 30.4, 15.5, 4.3};
		//管段直径(m)
		double[] dpl = new double[]{0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3};
		//摩|阻系数
		double[] slp = new double[]{0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017};
		// 上游管底高程(m)，下游管底高程(m)
		double[] ZJup = new double[]{3.894, 3.842, 3.784, 3.733, 3.687, 3.643, 3.601, 3.731, 3.886};
		double[] ZJdw = new double[]{3.842, 3.784, 3.733, 3.687, 3.643, 3.601, 3.546, 3.171, 3.7};
/**
		节点序号	1	2	3	4	5	6	7   8
		路径1	    1	2	3	4	5	6	7   8
		路径2	    9	8	-99
		路径3	    10	7	-99
*/
		// 管网路径数和路径节点号(-99表示空节点)
		int[][] Mroute = new int[][]{{1,2,3,4,5,6,7,8},{9,8,-99},{10,7,-99}};
		
/**
 *    	节点序号	1	2	3	4	5	6   7
		路径1	    7	6	5	4	3	2	1
		路径2	    8   -99
		路径3		9   -99

 */
		// 子系统分支路径管段数据矩阵 倒序pipe branches-reverse order
		int[][] Mbranch = new int[][]{{7,6,5,4,3,2,1},{8,-99},{9,-99}};


		
		/**
		 * 程序开始
		 * */		
		int ip, it, k1, kp, in1, in2, in3;
		int i, j, k, ik, jk, k00, kkk1 = 1;
		double dtnt, taa, tbb, AA, XX1, XX2;
		
		// XX[NT],qit[NT],sumqj[Nbr] [NT][NN],sumAj[Nbr] [NT][NN],Tnode[Nbr]
		// [NN][NN],sumTnode[Nbr] [NN][NN];
		double[] XX = new double[NT];			//时段降雨强度计算中间值 qit[it] = 167.0 * XX[it] / 1000.0;
		double[] qit = new double[NT];			//时段降雨强度（m3/min）
		double[][] sumqj = new double[NT][NN];	//时段节点受流量（m3/s）
		double[][] sumAj = new double[NT][NN];	//时段节点集水面积（m2）
		double[][] Tnode = new double[NN][NN];	//相邻节点间流经时间
		double[][] sumTnode = new double[NN][NN];//不相邻节点间总流经时间

		double hdc_min, hdc_max, fun_hd0, hda, hdb, hdc = 1.0, hde = 1.0, hdf, hdd;
		// hdcc[2],fun_hd[2],fun_hd0,hdcc0[NT][NP],hda,hdb,hdc,hde,hdf,hdd;
		double[] hdcc = new double[2];
		double[] fun_hd = new double[2];
		double[][] hdcc0 = new double[NT][NP];

		// 管段流速（m/s）
		// double vp[Nbr] [NP],slop[Nbr] [NP],Ad0,AD,RD,vp0=1.0;
		double[] vp = new double[NP];
		double[] slop = new double[NP];
		double Ad0, AD, RD, vp0 = 1.0;

		// double qpt[Nbr] [NT][NP],vpt[Nbr] [NT][NP],rid[Nbr]
		// [NT][NP],slopt[Nbr] [NT][NP],
		// Hwup[Nbr] [NT][NP],Hwdw[Nbr] [NT][NP],qjt[Nbr] [NT][NN],
		// overflow[Nbr] [NT][NN],Hw_over[Nbr]
		// [NT][NN],hdj0,ARD23,totalflow[Nbr] [NN],totalHw[Nbr] [NN];
		double[][] qpt = new double[NT][NP];//时段管段流量（m3/s）
		double[][] vpt = new double[NT][NP];//时段管段流速（m/s）
		double[][] rid = new double[NT][NP];//时段水力半径
		double[][] slopt = new double[NT][NP];//时段水力坡度
		double[][] Hwup = new double[NT][NP];//时段上端水位
		double[][] Hwdw = new double[NT][NP];//时段下端水位
		double[][] qjt = new double[NT][NN];//时段节点流量（m3/s）
		double[][] overflow = new double[NT][NN];//时段节点溢流量（m3/s）
		double[][] Hw_over = new double[NT][NN];//时段节点地面积水深度
		double hdj0, ARD23 = 1;
		double[] totalflow = new double[NN];//节点总积水量（m3）
		double[] totalHw = new double[NN];

		// 暴雨公式shanghai storm water formular:
		// (A1+C*lgP)/(t+b)**n---ln(N)=2.303log(N)---出口水位（m）
		double A1 = 17.53, C_storm = 0.95, b_storm = 11.77, P_simu = 10, n_storm = 0.88, dt = 2.0, rc = 0.375;
		double Hw_end = 1;



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
		
		// ================= 赋初值 ===============
		// for (kk = 0; kk < Nbr; kk++) {
		for (i = 0; i < NT; i++) {
			for (j = 0; j < NN; j++) {
				sumAj[i][j] = 0;
			}
		}
		for (i = 0; i < NT; i++) {
			for (j = 0; j < NN; j++) {
				sumqj[i][j] = 0;
			}
		}
		for (i = 0; i < NN; i++) {
			for (j = 0; j < NN; j++) {
				if (i == j) {
					Tnode[i][j] = 0;
				} else {
					Tnode[i][j] = -99;
				}
			}
		}
		for (i = 0; i < NN; i++) {
			for (j = 0; j < NN; j++) {
				sumTnode[i][j] = 0;
			}
		}
		// }
		// ==================Tnode-sumTnode=========================
		// for (kk = 0; kk < Nbr; kk++) {
		for (i = 0; i < NP; i++) {
			vp[i] = vp0;
		}
		for (kp = 0; kp < NP; kp++) {
			in1 = I0[kp]-1;
			in2 = J0[kp]-1;
			Tnode[in1][in2] = lp[kp] / vp[kp] / 10;
			slop[kp] = (ZJup[kp] - ZJdw[kp]) / lp[kp];
		}
		//
		for (i = 0; i < Nroute; i++) {
			for (j = 0; j < Nr_node; j++) {
				in1 = Mroute[i][j]-1;
				if (in1 >= 0) {
					for (k = j + 1; k < Nr_node; k++) {
						in2 = Mroute[i][k - 1]-1;
						in3 = Mroute[i][k]-1;
						if (in3 >= 0) {
							sumTnode[in1][in3] = sumTnode[in1][in2]
									+ Tnode[in2][in3];
						}
						else
							break;
					}
				}
				else
					break;
			}
		}
		// }
		// kk
		// =====print Mroute[i][j], Tnode, sumTnode,Mbranch[i][j]====
		// cout<<"pipe no.  I0    J0"<<endl;
		printStream.println("pipe no.  I0     J0");
		// for (kk = 0; kk < Nbr; kk++) {
		for (i = 0; i < NP; i++) {
			printStream.printf("%6d%6d%6d", i,I0 [i],J0 [i]);
			// s.Format("%6d%6d%6d\n", i,I0 [i],J0 [i]);
			// cout<<s;
			printStream.println();
		}
		// outfile<<endl;
		//
		// outfile<<" ip=";
		printStream.print("ip=");
		for (i = 0; i < NP; i++) {
			printStream.printf("%4d",i);
			// s.Format("%4d",i);
			// outfile<<s;
		}
		printStream.println();
		// outfile<<endl;
		//
		// outfile<<" I0=";
		printStream.print("I0=");
		for (i = 0; i < NP; i++) {
			printStream.printf("%4d",I0 [i]);
			// s.Format("%4d",I0 [i]);

			// outfile<<s;
		}
		printStream.println();
		// outfile<<endl;
		//
		// outfile<<" J0=";
		printStream.print("J0=");
		for (i = 0; i < NP; i++) {
			printStream.printf("%4d",J0 [i]);
			// s.Format("%4d",J0 [i]);
			// outfile<<s;
		}
		printStream.println();
		// outfile<<endl;
		// }// kk
		//
		// outfile<<endl;
		// outfile<<"===========  print Mroute [i][j]"<<endl;
		printStream.println("===========  print Mroute [i][j]");
		// for (kk = 0; kk < Nbr; kk++) {
		for (i = 0; i < Nroute; i++) {
			for (j = 0; j < Nr_node; j++) {
				printStream.printf("%6d",Mroute [i][j]);
				// s.Format("%6d",Mroute [i][j]);
				// outfile<<s;
			}
			printStream.println();
			// outfile<<endl; }
			// }// kk
			//
			// outfile<<endl;
			// outfile<<"===========  print Mbranch [i][j]"<<endl;
			printStream.println("===========  print Mbranch[i][j]:");
			// for (kk = 0; i < Nbr; kk++) {
			for (i = 0; i < Nstart; i++) {
				for (j = 0; j < Npline; j++) {
					printStream.printf("%6d",Mbranch [i][j]);
					// s.Format("%6d",Mbranch [i][j]);
				}
				printStream.println();
				// outfile<<endl; }
				// }// kk
				//
				// outfile<<"===========  print Tnode [i][j]"<<endl;
				// outfile<<"====j=  "<<endl;
				// outfile<<"      ";
				printStream.println("===========  print Tnode[i][j]:");
				printStream.print("==j=  ");
				// for (kk = 0; i < Nbr; kk++) {
				for (j = 0; j < NN; j++) { // s.Format("%6d",j);
											// outfile<<s;
					printStream.printf("%6d",j);
				}
				printStream.println();
				// outfile<<endl;
				//
				for (i = 0; i < NN; i++) {
					if (i < 10) {
						// outfile<<"i="<<i<<"   ";
						printStream.print("i=" + i + "   ");
					} else {
						// outfile<<"i="<<i<<"  ";
						printStream.print("i=" + i + "   ");

					}
					for (j = 0; j < NN; j++) {
						if (Tnode[i][j] < 0.0) {
							// outfile<<"      ";
							 printStream.print("      ");
						} else {
							// s.Format("%6.2f",Tnode [i][j]);
							printStream.printf("%6.2f",Tnode [i][j]);
							// outfile<<s;
						}
					}
					printStream.println();
					// outfile<<endl;
				}
				// }// kk
				//
				// outfile<<endl;
				// outfile<<"===========  print sumTnode [i][j]"<<endl;
				// outfile<<"==j=  ";
				// outfile<<"      ";
				printStream.println("\n===========  print sumTnode[i][j]:");
				printStream.print("==j=  ");
				// for (kk = 0; i < Nbr; kk++) {
				for (j = 0; j < NN; j++) { // s.Format("%6d",j);
											// outfile<<s;
					printStream.printf("%6d", j);
				}
				
				// outfile<<endl;
				//
				printStream.println();
				for (i = 0; i < NN; i++) {// outfile<<"i="<<i<<"   ";

					printStream.print("i=" + i + "   ");
					for (j = 0; j < NN; j++) {
						if (sumTnode[i][j] <= 0.0) {// outfile<<"      ";

							printStream.print("      ");
						} else

						{ // s.Format("%6.2f",sumTnode [i][j]);
							// outfile<<s;
							printStream.printf("%6.2f",sumTnode [i][j]);

						}
					}
					printStream.println();
					// outfile<<endl;
				}
				// }// kk
				// ================= 管网稳态流动模拟===================
				//
				// -----------动态模拟流量计算------------------------------------
				// ----------------节点汇水面积(ha)和汇水流量(m3/sec)计算--------
				// outfile<<endl;
				// outfile<<"===========  管网动态模拟计算      重现期＝ "<<P_simu<<"  年   时段数＝ "<<NT<<"       终点水位＝ "<<Hw_end<<"  m  ========="<<endl;
				printStream.println("\n===========  管网动态模拟计算      重现期＝ "
						+ P_simu + "  年   时段数＝ " + NT + "       终点水位＝ "
						+ Hw_end + "  m  =========");
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
				printStream.println("    it      dtnt      XX[it]     qit[it]");
				for (it = 0; it < NT; it++) {
					dtnt = dt * (float) (it);
					// s.Format("%6d%10.2lf%12.6lf%12.6lf\n",it,dtnt,XX[it],qit[it]);
					// outfile<<s;
					printStream.printf("%6d%10.2f%12.6f%12.6f\n",it,dtnt,XX[it],qit[it]);
					printStream.println();
				}
				// outfile<<endl;
				// xxxxxxx
				//
				// for (kk = 0; i < Nbr; kk++) {
				for (it = 0; it < NT; it++) {
					dtnt = dt + dt * (float) (it);
					for (j = 0; j < NN; j++) {
						sumAj[it][j] = Aj[j];
						sumqj[it][j] = Aj[j] * qit[it] * Acoef[j];
						for (i = 0; i < NN; i++) {
							if (sumTnode[i][j] > 0 && sumTnode[i][j] < dtnt) {
								sumAj[it][j] = sumAj[it][j] + Aj[i];
								sumqj[it][j] = sumqj[it][j] + Aj[i] * qit[it]
										* Acoef[i];
							}
						}
					}
				}
				// print sumAj[it][j] and sumqj[it][j]
				// outfile<<"  sumAj [it][j]="<<endl;
				printStream.println("  sumAj [it][j]=");
				for (it = 0; it < NT; it++) {
					for (j = 0; j < NN; j++) {// s.Format("%8.2lf",sumAj
												// [it][j]);
												// outfile<<s;
						printStream.printf("%8.2f",sumAj[it][j]);
					}
					printStream.println();
					// outfile<<endl;
				}
				printStream.println();
				// outfile<<endl;
				//
				// outfile<<"  sumqj [it][j]="<<endl;
				printStream.println("  sumqj [it][j]=");
				for (it = 0; it < NT; it++) {
					for (j = 0; j < NN; j++) {// s.Format("%8.2lf",sumqj
												// [it][j]);
												// outfile<<s;
						printStream.printf("%8.2f",sumqj[it][j]);
					}
					printStream.println();
					// outfile<<endl;
				}
				printStream.println();
				// outfile<<endl;
				// }// kk
				// ---------------------------------------------------------------
				for (it = 0; it < NT; it++) {
					for (i = 0; i < NN; i++) {
						overflow[it][i] = 0.0;
						Hw_over[it][i] = 0.0;
					}
				}
				for (it = 0; it < NT; it++) {
					for (j = 0; j < NP; j++)
						qpt[it][j] = -99.0;
				}
				// ---------------------------------------------------------------
				// kk-----000开始
				// for (kk = 0; i < Nbr; kk++) {
				for (it = 0; it < NT; it++)
				// --1
				{// outfile<<" it="<<it<<"  qpt [it][k]=";
					printStream.print(" it=" + it + "  qpt [it][k]=");
					for (j = 0; j < NN; j++) {
						for (k = 0; k < NP; k++) {
							if (I0[k]-1 == j) {
								qpt[it][k] = sumqj[it][j];
								// s.Format("%8.2lf",qpt [it][k]);
								// outfile<<s;
								printStream.printf("%8.2f",qpt [it][k]);
							}
						}
					}
					printStream.println();
					// outfile<<endl;
					// -------------------???????????????090127?????????????????------------------------
					for (ik = 0; ik < Nstart; ik++)
					// --2
					{
						for (jk = 0; jk < Npline; jk++)
						// --3
						{
							kp = Mbranch[ik][jk]-1;
							if (kp >= 0)
							// --4
							{
								if (J0[kp] == Nend) {
									Hwdw[it][kp] = Hw_end;
									// outfile<<"   it= "<<it<<"   kp= "<<kp<<"  Hwdm= "<<Hwdw[it][kp]<<"  Hw_end= "<<Hw_end<<endl;
									printStream.println("   it= " + it
											+ "   kp= " + (kp+1) + "  Hwdm= "
											+ df.format(Hwdw[it][kp]) + "  Hw_end= "
											+ Hw_end);
								} else {
									for (k1 = 0; k1 < NP; k1++) {
										if (I0[k1] == J0[kp])
											Hwdw[it][kp] = Hwup[it][k1];
									}
								}
								//
								Ad0 = 0.7854 * Math.pow(dpl[kp], 2.0);
								hdj0 = ZJdw[kp] + dpl[kp];
								if (Hwdw[it][kp] >= hdj0) {
									hdcc0[it][kp] = 1.0;
									rid[it][kp] = dpl[kp] / 4.0;
									vpt[it][kp] = qpt[it][kp] / Ad0;
									slopt[it][kp] = 10.29
											* Math.pow(slp[kp], 2.0)
											* Math.pow(qpt[it][kp], 2.0)
											/ Math.pow(dpl[kp], 5.333);
									Hwup[it][kp] = Hwdw[it][kp] + slopt[it][kp]
											* lp[kp];
									if (Hwup[it][kp] >= Hj[kp]) {
										Hwup[it][kp] = Hj[kp];
										slopt[it][kp] = (Hwup[it][kp] - Hwdw[it][kp])
												/ lp[kp];
										vpt[it][kp] = Math.pow(rid[it][kp],
												0.6667)
												* Math.pow(slopt[it][kp], 0.5)
												/ slp[kp];
										qpt[it][kp] = vpt[it][kp] * Ad0;
									}
								} else
								// --5
								{
									hdc_min = (Hwdw[it][kp] - ZJdw[kp])
											/ dpl[kp];
									slopt[it][kp] = slop[kp];
									if (hdc_min < 0.0)
										hdc_min = 0.0;
									ARD23 = slp[kp] * qpt[it][kp]
											/ Math.pow(slopt[it][kp], 0.5);
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
											hdc = Math.pow(Math.cos(hda), -1);
											hdd = hdcc[k00] * hdb;
											hde = hda * Math.pow(hdd, 0.5);
											hdf = Math.pow(dpl[kp], 2.0);
											AD = hdf / 4.0 * hdc - hdf / 2.00
													* hde;
											RD = dpl[kp] / 4.0 - dpl[kp] * hde
													/ (2.0 * hdc);
											fun_hd[k00] = AD
													* Math.pow(RD, 2.0 / 3.0)
													- ARD23;
										}
										if (Math.abs(fun_hd[0]) > Math.abs(fun_hd[1])) {
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
											printStream.println("it=  " + it
													+ "  ik=  " + ik + "jk=  "
													+ jk + "  kp=  " + (kp+1)
													+ "  kkk1=  " + kkk1
													+ "  hdcc[0]=  " + df.format(hdcc[0])
													+ "  hdcc[1]=  " + df.format(hdcc[1])
													+ "  fun_hd[0]= "
													+ df.format(fun_hd[0])
													+ "  fun_hd[1]= "
													+ df.format(fun_hd[1]) + "  ARD23= "
													+ ARD23);
										}
									}// 6-kkk1
										//
									rid[it][kp] = (0.25 - 0.5 * hde / hdc)
											* dpl[kp];
									Hwdw[it][kp] = ZJdw[kp] + hdcc0[it][kp]
											* dpl[kp];
									vpt[it][kp] = Math.pow(rid[it][kp],
											2.0 / 3.0)
											* Math.pow(slopt[it][kp], 0.5)
											/ slp[kp];
									Hwup[it][kp] = Hwdw[it][kp] + slopt[it][kp]
											* lp[kp];
									// ---------h/d---0.618aaa---------------------------
								} // 5--end_hdcc0[it][kp]
									//
									// outfile<<"   it= "<<it<<"   kp= "<<kp<<"   I0[kp]= "<<I0
									// [kp]<<"  Hwdm= "<<Hwdw[it][kp]<<"  Hj= "<<Hj[I0
									// [kp]]
									// <<"  hdcc0= "<<hdcc0[it][kp]<<"  qpt= "<<qpt[it][kp]<<"  vpt= "<<vpt[it][kp]<<"  hdcc0= "<<hdcc0[it][kp]<<" fun_hd0= "<<fun_hd0<<endl;
								
								printStream.println("   it= " + it + "   kp= " + kp 
										+ "   I0[kp]= " + df.format(I0[kp]) + "  Hwdm= " + df.format(Hwdw[it][kp]) 
										+ "  Hj= " + df.format(Hj[I0[kp]-1]) + "  hdcc0= " + df.format(hdcc0[it][kp]) 
										+ "  qpt= " + df.format(qpt[it][kp]) + "  vpt= " + df.format(vpt[it][kp]) 
										+ "  hdcc0= " + df.format(hdcc0[it][kp])+ " fun_hd0= "+ df.format(fun_hd[0]));

							}// 4 if(kp>=0) end
							else
								break;
						}// 3 ---jk end ---
					}// --2---ik end ---Hj[kp]
						// --------------------------------- 开始计算溢流节点
						// ---------------
						//
					for (i = 0; i < NP; i++) 
					{
						j = I0[i]-1;
						//cj20160828 疑问  时段上端水位大于节点地面标高？
						if (Hwup[it][i] == Hj[j])
						{
							overflow[it][j] = overflow[it][j] - qpt[it][i];
							for (ip = 0; ip < NP; ip++) 
							{
								k1 = J0[ip]-1;
								if (k1 == j)
									overflow[it][j] = overflow[it][j] + qpt[it][ip];
							}
							qjt[it][j] = Aj[j] * qit[it] * Acoef[j];
							overflow[it][j] = (overflow[it][j] + qjt[it][j]) * dt * 60.0;
							Hw_over[it][j] = overflow[it][j] / Aj[j] / 10000.0 * 1000.0;
						}
					}
					// ----qjt[it][j],overflow[it][j],totalflow[j],totalHw[j]
					// ------------------ 计算溢流节点结束 ---------------
					// outfile<<endl;
					// outfile<<"    it   管段号  I0   J0   管径dpl   管段qp 水力半径R    充满度¨¨ 流速¨′(m/s)  上游水位  下游水位  上管底高  下管底高  上地面高"<<endl;
					printStream.println("    it   管段号  I0   J0   管径dpl   管段qp 水力半径R    充满度 流速(m/s)  上游水位  下游水位  上管底高  下管底高  上地面高");
					for (i = 0; i < NP; i++) {// s.Format("%6d%6d%6d%5d%10.2lf%10.3lf%10.3lf%10.3lf%10.3lf%10.3lf%10.3lf%10.3lf%10.3lf%10.3lf\n",
												// it,i,I0 [i],J0
												// [i],dpl [i],qpt
												// [it][i],rid
												// [it][i],hdcc0
												// [it][i],vpt
												// [it][i],Hwup
												// [it][i],Hwdw
												// [it][i],ZJup
												// [i],ZJdw [i],Hj
												// [I0[i]]);
						// outfile<<s;
						printStream.printf("%6d%6d%6d%5d%10.2f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f",
								 it,i,I0 [i],J0[i],dpl [i],qpt[it][i],rid[it][i],hdcc0[it][i],
								 vpt[it][i],Hwup[it][i],Hwdw[it][i],ZJup[i],ZJdw [i],Hj[I0[i]-1]);
						printStream.println();
					}
					// outfile<<endl;
					printStream.println();
				}// 1-- it end ---
					// ------- 计算节点积水量和积水深度¨¨(m)----
				for (j = 0; j < NN; j++) {
					totalflow[j] = 0.0;
					totalHw[j] = 0.0;
				}
				for (j = 0; j < NN; j++) {
					for (it = 0; it < NT; it++) {
						totalflow[j] = totalflow[j] + overflow[it][j];
						totalHw[j] = totalHw[j] + Hw_over[it][j];
					}
				}
				// -----屏幕输出管段水力计算结束------
				// cout<<" it 管段号 I0 J0 管径dpl 管段qp 水力半径R 充满度 流速(m/s) 上游水位
				// 下游水位 上管底高 下管底高 上地面高”<<endl;
				printStream.println("    it   管段号  I0   J0   管径dpl   管段qp 水力半径R    充满度 流速(m/s)  上游水位  下游水位  上管底高  下管底高  上地面高");
				// outfile<<endl;
				for (it = 0; it < NT; it++) {// outfile<<"  it= "<<it<<endl;
												// outfile<<endl;
												// outfile<<" it 管段号 I0 J0
												// 管径dpl 管段qp 水力半径R 充满度¨¨
												// 流速(m/s) 上游水位 下游水位 上管底高程
												// 下管底高程<<endl;
					for (i = 0; i < NP; i++) {// s.Format("%6d%6d%6d%5d%10.2lf%10.3lf%10.3lf%10.3lf%10.3lf%10.3lf%10.3lf%10.3lf%10.3lf%10.3lf\n",
												// it,i,I0 [i],J0
												// [i],dpl [i],qpt
												// [it][i],rid
												// [it][i],hdcc0
												// [it][i],vpt
												// [it][i],Hwup
												// [it][i],Hwdw
												// [it][i],ZJup
												// [i],ZJdw [i],Hj
												// [I0[i]]);
												// outfile<<s;
						// cout<<s;
						printStream.printf("%6d%6d%6d%5d%10.2f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f",
								 it,i,I0 [i],J0[i],dpl [i],qpt[it][i],rid[it][i],
								 hdcc0[it][i],vpt[it][i],Hwup[it][i],Hwdw[it][i],ZJup[i],ZJdw [i],Hj[I0[i]-1]);
						printStream.println();
					}
					// outfile<<endl;
					printStream.println();
				}
				// outfile<<endl;
				// ------------- 节点溢流计算结果输出 ---------------
				// outfile<<" ======== 时段节点积水量(m3) ========"<<endl;
				// outfile<<"  i=    ";
				printStream.println(" ======== 时段节点积水量(m3) ========");
				printStream.print("  i=    ");
				for (i = 0; i < NN; i++) {
					if (i < 10) {
						// outfile<<"  "<<i<<"   ";
						printStream.print(" " + i + "   ");
					} else {
						// outfile<<" "<<i<<"   ";
						printStream.print(i + "   ");
					}
				}				
				// outfile<<endl;
				// outfile<<"it="<<endl;
				printStream.println();
				printStream.println("it=");
				for (it = 0; it < NT; it++) {
					if (it < 10) {
						// outfile<<" "<<it<<"   ";
						printStream.print(" " + it + "   ");
					} else {
						// outfile<<it<<"   ";
						printStream.print(it + "   ");
					}
					//
					for (i = 0; i < NN; i++) {
						if (overflow[it][i] <= 0.0) {
							// outfile<<"      ";
							printStream.print("      ");
						} else {
							// s.Format("%6.1f",overflow [it][i]);
							// outfile<<s;
							printStream.printf("%6.1f",overflow [it][i]);
						}
					}
					// outfile<<endl;
					printStream.println();
				}
				//
				// outfile<<" ======== 时段节点积水深度输出(mm) ========"<<endl;
				// outfile<<"  i=    ";
				printStream.println(" ======== 时段节点积水深度输出(mm) ========");
				printStream.print("  i=    ");
				for (i = 0; i < NN; i++) {
					if (i < 10) {
						// outfile<<"  "<<i<<"   ";
						printStream.print("  " + i + "  ");
					} else {
						// outfile<<" "<<i<<"   ";
						printStream.print(i + "  ");
					}
				} // outfile<<endl;
					// outfile<<"it="<<endl;
				printStream.println();
				printStream.println("it=");
				for (it = 0; it < NT; it++) {
					if (it < 10) {
						// outfile<<" "<<it<<"   ";
						printStream.print(" " + it + "  ");
					} else {
						// outfile<<it<<"   ";
						printStream.print(it + "  ");
					}
					//
					for (i = 0; i < NN; i++) {
						if (Hw_over[it][i] <= 0.0) {
							// outfile<<"      ";
							printStream.print("      ");
						} else {
							// s.Format("%6.1f",Hw_over [it][i]);
							// outfile<<s;
							printStream.printf("%6.1f",Hw_over [it][i]);
						}
					}
					printStream.println();
					// outfile<<endl;
				}
			}// kk---000结束
				//
				// outfile.close();

		}
		printStream.close();
		// }
	}
}
