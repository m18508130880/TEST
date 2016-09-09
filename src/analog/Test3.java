package analog;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.text.DecimalFormat;

//20150908new��ˮ��������ģ��Ƽ������-֥�Ӹ������-�ٽ�ˮ���㷨
//����SQ LIU, TONGJI UNIVERSITY, 8 SEP 2016 
//
public class Test3 {
	public static void main(String[] args) throws FileNotFoundException {

		String FileName = "201608031-��-��ˮ��������ģ��-֥�Ӹ������-���ҳ�002.txt";
		FileOutputStream fs = new FileOutputStream(new File(FileName));
		PrintStream printStream = new PrintStream(fs);
		printStream.println(FileName);

		DecimalFormat df = new DecimalFormat("##.####");

		// CString s;
		// �����������ݣ�
		// �ܶ������ڵ������ܵ��������·�����ܶ����������������ģ��ʱ������֥�Ӹ���ʱ��λ��
		// �ܵ�·������·�����ڵ������յ�ڵ�ţ��м�������ļ�ָ��
		// const int
		// NP=9,NN=10,Nstart=3,Npline=7,NT=60,NR=23,Nroute=3,Nr_node=8,Nend=7,Iprt=0;
		int NP = 9, NN = 10, Nstart = 3, Npline = 7;
		int NT = 60;
		int NR = 23, Nroute = 3, Nr_node = 8, Nend = 7, Iprt = 0;
		// ���깫ʽ����shanghai storm water formular:
		// (A1+C*lgP)/(t+b)**n---ln(N)=2.303log(N)---����ˮλ��m��
		// �ܶ����٣�m/s��, �ܶ��趨����vp0�����氼͹ϵ��csf
		double A1 = 17.53, C_storm = 0.95, tmin = 10, b_storm = 11.77, P_simu = 50, n_storm = 0.88, dt = 2.0, rc = 0.375, Hw_end = 3.0, vp0 = 0.8, csf = 3.0;
		// �ڵ��ˮ���(ha)
		double[] Aj = new double[] { 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2,
				0.2, 0.2 };
		// �ڵ��ˮ������ϵ��
		double[] Acoef = new double[] { 0.62, 0.62, 0.62, 0.62, 0.62, 0.62,
				0.62, 0.62, 0.62, 0.62 };
		// �ڵ�����ߣ�m��
		double[] Hj = new double[] { 5.24, 5.19, 5.18, 5.20, 5.21, 5.20, 5.20,
				5.12, 5.13, 5.18 };
		// ����·������·���ڵ��(��99��ʾ�սڵ�)
		int[][] Mroute = new int[][] { { 0, 1, 2, 3, 4, 5, 6, 7 },
				{ 8, 7, -99, -99, -99, -99, -99, -99 },
				{ 9, 6, -99, -99, -99, -99, -99, -99 } };

		// pipe branches-reverse order
		int[][] Mbranch = new int[][] { { 6, 5, 4, 3, 2, 1, 0 },
				{ 7, -99, -99, -99, -99, -99, -99, },
				{ 8, -99, -99, -99, -99, -99, -99 } };
		// �ܶ����νڵ��I0,���νڵ��J0���ܶγ���(m),Ħ��ϵ��
		int[] I0 = new int[] { 0, 1, 2, 3, 4, 5, 6, 8, 9 };
		int[] J0 = new int[] { 1, 2, 3, 4, 5, 6, 7, 7, 6 };
		double[] lp = new double[] { 50, 50, 50, 50, 50, 50, 50, 50, 50 };
		double[] slp = new double[] { 0.014, 0.014, 0.014, 0.014, 0.014, 0.014,
				0.014, 0.014, 0.014 };
		// double
		// slp[NP]={0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013,0.013};
		// �ܶ�ֱ��(m)�����ιܵ׸߳�(m)�����ιܵ׸߳�(m)
		double[] dpl = new double[] { 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3,
				0.3 };
		double[] ZJup = new double[] { 3.89, 3.84, 3.78, 3.73, 3.68, 3.64,
				3.60, 3.73, 3.88 };
		double[] ZJdw = new double[] { 3.84, 3.78, 3.73, 3.68, 3.64, 3.60,
				3.55, 3.60, 3.70 };
		// ----�ٽ�ˮ��������----
		double sita0 = 3.0, eps = 0.001, alfa = 0.5;
		double Ad0, qkpmax, sita, cons_b, sita_s = 0, sita_c, fsita, dfdsita, dfsita, ssita = 0, csita = 0, hyd_A, hafsita, shafsita = 0, chafsita, sita_p = 0;
		// �м����
		int i, j, k = 0, ik, jk, it, k1, kp, in1, in2, in3;
		double dtnt, taa, tbb, AA, XX1, XX2, hdj0;
		double[] XX = new double[NT];
		double[] qit = new double[NT];
		double[][] sumqj = new double[NT][NN];
		double[][] sumAj = new double[NT][NN];
		double[][] Tnode = new double[NN][NN];
		double[][] sumTnode = new double[NN][NN];
		double[] vp = new double[NP];
		double[] slop = new double[NP];
		double[][] qpt = new double[NT][NP];
		double[][] qqkp = new double[NT][NP];
		double[][] vpt = new double[NT][NP];
		double[][] rid = new double[NT][NP];
		double[][] slopt = new double[NT][NP];
		double[][] Hwup = new double[NT][NP];
		double[][] Hwdw = new double[NT][NP];
		double[][] hdcc0 = new double[NT][NP];
		double[][] overflow = new double[NT][NN];
		double[][] Hw_over = new double[NT][NN];
		//
		// ----------------------------------------------------------------------------------------------------------
		//
		// --��������ļ���ʼ---
		// ofstream outfile;
		// outfile.open("20160908-��ˮ��������ģ��-���ҳ�-50-3.txt");
		printStream.print("20160908-��ˮ��������ģ��-���ҳ�-50-3.txt");
		//
		// cout<<"------ ��ˮ��������ģ��-���ҳ� ------"<<endl;
		System.out.println("------ ��ˮ��������ģ��-���ҳ� ------");
		// ================= ����ֵ ===============================
		for (i = 0; i < NT; i++) {
			for (j = 0; j < NN; j++)
				sumAj[i][j] = 0;
		}
		for (i = 0; i < NT; i++) {
			for (j = 0; j < NN; j++)
				sumqj[i][j] = 0;
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
			for (j = 0; j < NN; j++)
				sumTnode[i][j] = 0;
		}
		// ==================Tnode-sumTnode=========================
		for (i = 0; i < NP; i++)
			vp[i] = vp0;
		for (kp = 0; kp < NP; kp++) {
			in1 = I0[kp];
			in2 = J0[kp];
			Tnode[in1][in2] = lp[kp] / vp[kp] / 60;
			slop[kp] = (ZJup[kp] - ZJdw[kp]) / lp[kp];
		}
		//
		for (i = 0; i < Nroute; i++) {
			for (j = 0; j < Nr_node; j++) {
				in1 = Mroute[i][j];
				if (in1 >= 0) {
					for (k = j + 1; k < Nr_node; k++) {
						in2 = Mroute[i][k - 1];
						in3 = Mroute[i][k];
						if (in3 >= 0) {
							sumTnode[in1][in3] = sumTnode[in1][in2]
									+ Tnode[in2][in3];
						}
					}
				}
			}
		}
		// =====print Mroute[i][j], Tnode, sumTnode,Mbranch[i][j]====
		// cout<<"pipe no.  I0    J0"<<endl;
		System.out.println("pipe no.  I0    J0");
		for (i = 0; i < NP; i++) {// s.Format("%6d%6d%6d\n", i,I0[i],J0[i]);
									// cout<<s;
			System.out.printf("%6d%6d%6d", i, I0[i], J0[i]);
			System.out.println();
		}
		// outfile<<endl;
		printStream.println();

		// outfile<<" ip=";
		printStream.print(" ip=");
		for (i = 0; i < NP; i++) {// s.Format("%4d",i);
									// outfile<<s;
			printStream.printf("%4d", i);
		}
		// outfile<<endl;
		printStream.println();
		// outfile<<" I0=";
		printStream.print(" I0=");
		for (i = 0; i < NP; i++) { // s.Format("%4d",I0[i]);
									// outfile<<s;
			printStream.printf("%4d", I0[i]);
		}
		// outfile<<endl;
		printStream.println();
		// outfile<<" J0=";
		printStream.print(" J0=");
		for (i = 0; i < NP; i++) { // s.Format("%4d",J0[i]);
									// outfile<<s;
			printStream.printf("%4d", J0[i]);
		}
		// outfile<<endl;
		printStream.println();
		// outfile<<endl;
		printStream.println();
		// outfile<<"===========  print Mroute[i][j]"<<endl;
		printStream.println("===========  print Mroute[i][j]");
		for (i = 0; i < Nroute; i++) {
			for (j = 0; j < Nr_node; j++) { // s.Format("%6d",Mroute[i][j]);
											// outfile<<s;
				printStream.printf("%6d", Mroute[i][j]);
			}
			// outfile<<endl;
			printStream.println();
		}
		// outfile<<endl;
		printStream.println();
		// outfile<<"===========  print Mbranch[i][j]"<<endl;
		printStream.println("===========  print Mbranch[i][j]");
		for (i = 0; i < Nstart; i++) {
			for (j = 0; j < Npline; j++) { // s.Format("%6d",Mbranch[i][j]);
											// outfile<<s;
				printStream.printf("%6d", Mbranch[i][j]);
			}
			// outfile<<endl;
			printStream.println();
		}
		// outfile<<"===========  print Tnode[i][j]"<<endl;
		printStream.println("===========  print Tnode[i][j]");
		// outfile<<"====j=  "<<endl;
		printStream.println("====j=  ");
		// outfile<<"      ";
		printStream.println("      ");
		for (j = 0; j < NN; j++) { // s.Format("%6d",j);
									// outfile<<s;
			printStream.printf("%6d", j);
		}
		// outfile<<endl;
		printStream.println();
		for (i = 0; i < NN; i++) {
			if (i < 10) {// outfile<<"i="<<i<<"   ";
				printStream.print("i=" + i + "   ");
			} else {// outfile<<"i="<<i<<"  ";
				printStream.print("i=" + i + "  ");
			}
			for (j = 0; j < NN; j++) {
				if (Tnode[i][j] < 0.0) {// outfile<<"      ";
					printStream.print("      ");
				} else { // s.Format("%6.2f",Tnode[i][j]);
							// outfile<<s;
					printStream.printf("%6.2f", Tnode[i][j]);
				}
			}
			// outfile<<endl;
			printStream.println();
		}
		// outfile<<endl;
		printStream.println();
		// outfile<<"===========  print sumTnode[i][j]"<<endl;
		printStream.println("===========  print sumTnode[i][j]");
		// outfile<<"==j=  ";
		printStream.println("==j=  ");
		// outfile
		for (j = 0; j < NN; j++) { // s.Format("%6d",j);
									// outfile<<s;
			printStream.printf("%6d", j);
		}
		// outfile<<endl;
		printStream.println();
		for (i = 0; i < NN; i++) {// outfile<<"i="<<i<<"   ";
			printStream.print("i=" + i + "   ");
			for (j = 0; j < NN; j++) {
				if (sumTnode[i][j] <= 0.0) {// outfile<<"      ";
					printStream.print("      ");
				} else { // s.Format("%6.2f",sumTnode[i][j]);
							// outfile<<s;
					printStream.printf("%6.2f", sumTnode[i][j]);
				}
			}
			// outfile<<endl;
			printStream.println();
		}
		// ================= ����׼��̬����ģ��============================
		//
		// -------------------��̬ģ����������-----------------------------
		// ----------------�ڵ��ˮ���(ha)�ͻ�ˮ����(m3/sec)����--------
		// outfile<<endl;
		printStream.println();
		// outfile<<"===========  ������̬ģ�����      �����ڣ� "<<P_simu<<"  ��   ʱ������ "<<NT<<"       �յ�ˮλ�� "<<Hw_end<<"  m  ========="<<endl;
		printStream.println("===========  ������̬ģ�����      �����ڣ� " + P_simu
				+ "  ��   ʱ������ " + NT + "       �յ�ˮλ�� " + Hw_end
				+ "  m  =========");
		// xxxxxxx
		// ֥�Ӹ������--rainfall intensity at every time step--
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
				XX1 = AA * ((1.0 - n_storm) * taa / (1.0 - rc) + b_storm);
				XX2 = Math.pow((taa / (1.0 - rc) + b_storm), (n_storm + 1.0));
			}
			XX[it] = XX1 / XX2;
			qit[it] = 167.0 * XX[it] / 1000.0;
		}
		//
		// outfile<<endl;
		printStream.println();
		// outfile<<"    it      dtnt      XX[it]     qit[it]"<<endl;
		printStream.println("    it      dtnt      XX[it]     qit[it]");
		for (it = 0; it < NT; it++) {
			dtnt = dt * (float) (it);
			// s.Format("%6d%10.2lf%12.6lf%12.6lf\n",it,dtnt,XX[it],qit[it]);
			// outfile<<s;
			printStream.printf("%6d%10.2f%12.6f%12.6f", it, dtnt, XX[it],
					qit[it]);
			printStream.println();
		}
		// outfile<<endl;
		printStream.println();
		// xxxxxxx
		//
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
		// outfile<<"  sumAj[it][j]="<<endl;
		printStream.println("  sumAj[it][j]=");
		for (it = 0; it < NT; it++) {
			for (j = 0; j < NN; j++) {// s.Format("%8.2lf",sumAj[it][j]);
										// outfile<<s;
				printStream.printf("%8.2f", sumAj[it][j]);
			}
			// outfile<<endl;
			printStream.println();
		}
		// outfile<<endl;
		printStream.println();
		// outfile<<"  sumqj[it][j]="<<endl;
		printStream.println("  sumqj[it][j]=");
		for (it = 0; it < NT; it++) {
			for (j = 0; j < NN; j++) {// s.Format("%8.2lf",sumqj[it][j]);
										// outfile<<s;
				printStream.printf("%8.2f", sumqj[it][j]);
			}
			// outfile<<endl;
			printStream.println();
		}
		// outfile<<endl;
		printStream.println();
		// ---------------------------------------------------------------
		for (it = 0; it < NT; it++) {
			for (i = 0; i < NN; i++) {
				overflow[it][i] = 0.0;
				Hw_over[it][i] = 0.0;
			}
		}
		for (it = 0; it < NT; it++) {
			for (j = 0; j < NP; j++) {
				qpt[it][j] = -99.0;
				qqkp[it][j] = 0.0;
			}
		}
		// ---------------------------------------------------------------
		for (it = 0; it < NT; it++)
		// --1
		{// outfile<<" it="<<it<<"  qpt[it][k]=";
			printStream.print(" it=" + it + "  qpt[it][k]=");
			for (j = 0; j < NN; j++) {
				for (k = 0; k < NP; k++) {
					if (I0[k] == j) {
						qpt[it][k] = sumqj[it][j];
						// s.Format("%8.2lf",qpt[it][k]); outfile<<s;
						printStream.printf("%8.2f", qpt[it][k]);
					}
				}
			}
			// outfile<<endl;
			printStream.println();
			// -------------------20090127-sql����------------------------
			for (ik = 0; ik < Nstart; ik++)
			// --2
			{
				for (jk = 0; jk < Npline; jk++)
				// --3
				{
					kp = Mbranch[ik][jk];
					if (kp >= 0)
					// --4
					{
						if (J0[kp] == Nend) {
							Hwdw[it][kp] = Hw_end;
							// outfile<<"   it= "<<it<<"   kp= "<<kp<<"  Hwdm= "<<Hwdw[it][kp]<<"  Hw_end= "<<Hw_end<<endl;
							if(1 == Iprt){
								printStream.println("   it= " + it + "   kp= " + kp
										+ "  Hwdm= " + Hwdw[it][kp] + "  Hw_end= "
										+ Hw_end);
							}
						} else {
							for (k1 = 0; k1 < NP; k1++) {
								if (I0[k1] == J0[kp])
									Hwdw[it][kp] = Hwup[it][k1];
							}
						}
						//
						Ad0 = 0.7854 * Math.pow(dpl[kp], 2.0);
						hdj0 = ZJdw[kp] + dpl[kp];
						if (Hwdw[it][kp] >= hdj0) {// outfile<<"   it= "<<it<<"   kp= "<<kp<<"  Hwdm= "<<Hwdw[it][kp]<<"  ��û���� "<<endl;
							if(1 == Iprt){
								printStream.println("   it= " + it + "   kp= " + kp
										+ "  Hwdm= " + df.format(Hwdw[it][kp]) + "  ��û���� ");
							}
							hdcc0[it][kp] = 1.0;
							rid[it][kp] = dpl[kp] / 4.0;
							vpt[it][kp] = qpt[it][kp] / Ad0;
							slopt[it][kp] = 10.29 * Math.pow(slp[kp], 2.0)
									* Math.pow(qpt[it][kp], 2.0)
									/ Math.pow(dpl[kp], 5.333);
							Hwup[it][kp] = Hwdw[it][kp] + slopt[it][kp]
									* lp[kp];
							if (Hwup[it][kp] >= Hj[I0[kp]]) {
								Hwup[it][kp] = Hj[I0[kp]];
								slopt[it][kp] = (Hwup[it][kp] - Hwdw[it][kp])
										/ lp[kp];
								if (slopt[it][kp] < 0.0) {
									slopt[it][kp] = Math.abs(slopt[it][kp]);
								}
								vpt[it][kp] = Math.pow(rid[it][kp], 0.6667)
										* Math.pow(slopt[it][kp], 0.5)
										/ slp[kp];
								qqkp[it][kp] = vpt[it][kp] * Ad0;
								if (qqkp[it][kp] < 0.0) {
									qqkp[it][kp] = Math.abs(qqkp[it][kp]);
								}
							}
						} else
						// --5
						{// outfile<<"   it= "<<it<<"   kp= "<<kp<<"  Hwdw= "<<Hwdw[it][kp]<<"  ����û���� "<<endl;
							if(1 == Iprt){
								printStream.println("   it= " + it + "   kp= " + kp
										+ "  Hwdw= " + df.format(Hwdw[it][kp]) + "  ����û���� ");
							}
							// --20160907�޸Ŀ�ʼ---�����ٽ�ˮ���㷨-----------------------
							//
							qkpmax = 2.46 * Math.pow(dpl[kp], 2.5);
							if (qpt[it][kp] > qkpmax * 0.95) {// outfile<<"   it= "<<it<<"   kp= "<<kp<<"  qkpmax= "<<qkpmax<<"  ����û���ܳ��� "<<endl;
								if(1 == Iprt){
									printStream.println("   it= " + it + "   kp= "
											+ kp + "  qkpmax= " + qkpmax
											+ "  ����û���ܳ��� ");
								}								
								Hwdw[it][kp] = ZJdw[kp] + dpl[kp] * 1.1;
								// {
								hdcc0[it][kp] = 1.0;
								rid[it][kp] = dpl[kp] / 4.0;
								vpt[it][kp] = qpt[it][kp] / Ad0;
								slopt[it][kp] = 10.29 * Math.pow(slp[kp], 2.0)
										* Math.pow(qpt[it][kp], 2.0)
										/ Math.pow(dpl[kp], 5.333);
								Hwup[it][kp] = Hwdw[it][kp] + slopt[it][kp]
										* lp[kp];
								if (Hwup[it][kp] >= Hj[I0[kp]]) {
									Hwup[it][kp] = Hj[I0[kp]];
									slopt[it][kp] = (Hwup[it][kp] - Hwdw[it][kp])
											/ lp[kp];
									if (slopt[it][kp] < 0.0) {
										slopt[it][kp] = Math.abs(slopt[it][kp]);
									}
									vpt[it][kp] = Math.pow(rid[it][kp], 0.6667)
											* Math.pow(slopt[it][kp], 0.5)
											/ slp[kp];
									qqkp[it][kp] = vpt[it][kp] * Ad0;
									if (qqkp[it][kp] < 0.0) {
										qqkp[it][kp] = Math.abs(qqkp[it][kp]);
									}
								}
							} else {// outfile<<"   it= "<<it<<"   kp= "<<kp<<"  Hwdm= "<<Hwdw[it][kp]<<"  ����û�����ܳ��� "<<endl;
								if(1 == Iprt){
									printStream.println("   it= " + it + "   kp= "
											+ kp + "  Hwdm= " + df.format(Hwdw[it][kp])
											+ "  ����û�����ܳ��� ");
								}
								i = 0;
								sita = sita0;
								cons_b = 0.276843 * Math.pow(dpl[kp], 2.5)
										/ qpt[it][kp];
								if (Iprt == 1) {// outfile<<"   k= "<<k<<"   qpt[it][kp]= "<<qpt[it][kp]<<"   cons_b= "<<cons_b<<endl;
									printStream.println("   k= " + k
											+ "   qpt[it][kp]= " + qpt[it][kp]
											+ "   cons_b= " + cons_b);
								}
								/*
								 * L200: ssita=Math.sin(sita);
								 * csita=Math.cos(sita); hafsita= sita/2.0;
								 * shafsita=Math.sin(hafsita);
								 * chafsita=Math.cos(hafsita);
								 * sita_s=sita-Math.sin(sita);
								 * sita_c=1-Math.cos(sita); sita_p=Math.pow((1.0
								 * - chafsita),-0.5);
								 */
								// ==================================
								fsita = cons_b * sita_s - sita_p;
								dfsita = Math.abs(fsita);
								if (dfsita < eps) {
									hdcc0[it][kp] = (1 - Math.cos(sita / 2)) / 2;
									rid[it][kp] = 0.25 * dpl[kp]
											* (sita - Math.sin(sita)) / sita;
									vpt[it][kp] = Math.pow(rid[it][kp], 0.6667)
											* Math.pow(slop[kp], 0.5) / slp[kp];
								} else {
									dfdsita = cons_b * (1.0 - csita) + 0.25
											* Math.pow(sita_p, -1.0) * shafsita;
									sita = sita - alfa * fsita / dfdsita;
									if (Iprt == 1) {// outfile<<"   i= "<<i<<"   sita= "<<sita<<"   ssita= "<<ssita<<"   csita= "<<csita<<"   fsita= "<<fsita<<"   dfdsita= "<<dfdsita<<endl;
										printStream.println("   i= " + i
												+ "   sita= " + sita
												+ "   ssita= " + ssita
												+ "   csita= " + csita
												+ "   fsita= " + fsita
												+ "   dfdsita= " + dfdsita);
									}
									i = i + 1;
									ssita = Math.sin(sita);
									csita = Math.cos(sita);
									hafsita = sita / 2.0;
									shafsita = Math.sin(hafsita);
									chafsita = Math.cos(hafsita);
									sita_s = sita - Math.sin(sita);
									sita_c = 1 - Math.cos(sita);
									sita_p = Math.pow((1.0 - chafsita), -0.5);
								}
								//
								/*
								 * L100: hdcc0[it][kp]=(1-cos(sita/2))/2;
								 * rid[it]
								 * [kp]=0.25*dpl[kp]*(sita-sin(sita))/sita;
								 * vpt[it
								 * ][kp]=pow(rid[it][kp],0.6667)*pow(slop[kp
								 * ],0.5)/slp[kp];
								 */
							}// ---for(k=0;k<N;k++)����
								//
								// --20160907�޸Ľ���---�����ٽ�ˮ���㷨-----------------------
								//
							Hwdw[it][kp] = ZJdw[kp] + hdcc0[it][kp] * dpl[kp]
									* 1.1;
							Hwup[it][kp] = Hwdw[it][kp] + slop[kp] * lp[kp];
						} // 5--end
							// ------- ���it������ ----------
						if (Iprt == 1) {/*
										 * outfile<<"   it= "<<it<<"   kp= "<<kp<<
										 * "   I0[kp]= "
										 * <<I0[kp]<<"  Hwdm= "<<Hwdw
										 * [it][kp]<<"  Hwup= "
										 * <<Hwup[it][kp]<<"  Hj= "<<Hj[I0[kp]]
										 * <<
										 * "  hdcc0= "<<hdcc0[it][kp]<<"  qpt= "
										 * <<
										 * qpt[it][kp]<<"  qqkp= "<<qqkp[it][kp
										 * ]<<"  vpt= "<<vpt[it][kp]<<endl;
										 */
							printStream.println("   it= " + it + "   kp= " + kp
									+ "   I0[kp]= " + I0[kp] + "  Hwdm= "
									+ Hwdw[it][kp] + "  Hwup= " + Hwup[it][kp]
									+ "  Hj= " + Hj[I0[kp]] + "  hdcc0= "
									+ hdcc0[it][kp] + "  qpt= " + qpt[it][kp]
									+ "  qqkp= " + qqkp[it][kp] + "  vpt= "
									+ vpt[it][kp]);
						}
					}// --4 if(kp>=0) end
				}// --3 ---jk end
			}// --2---ik end
				//
			// outfile<<endl;
			printStream.println();
			// outfile<<"    it   �ܶκ�  I0   J0 �ܾ�dpl     �ܶ�qp   ˮ���뾶R  ������ ����(m/s)  ����ˮλ  ����ˮλ  �Ϲܵ׸�  �¹ܵ׸�  �ܶ��¶�  �ϵ����"<<endl;
			printStream
					.println("    it   �ܶκ�  I0   J0 �ܾ�dpl     �ܶ�qp   ˮ���뾶R  ������ ����(m/s)  ����ˮλ  ����ˮλ  �Ϲܵ׸�  �¹ܵ׸�  �ܶ��¶�  �ϵ����");
			for (i = 0; i < NP; i++) {/*
									 * s.Format(
									 * "%6d%6d%6d%5d%8.2lf%12.3lf%8.3lf%8.3lf%10.3lf%10.3lf%10.3lf%10.3lf%10.3lf%10.5lf%10.3lf\n"
									 * ,
									 * it,i,I0[i],J0[i],dpl[i],qpt[it][i],rid[it
									 * ]
									 * [i],hdcc0[it][i],vpt[it][i],Hwup[it][i],
									 * Hwdw [
									 * it][i],ZJup[i],ZJdw[i],slop[i],Hj[I0[i]])
									 * ; outfile<<s;
									 */
				printStream
						.printf("%6d%6d%6d%5d%8.2f%12.3f%10.3f%8.3f%10.3f%10.3f%10.3f%10.3f%10.3f%10.5f%10.3f",
								it, i, I0[i], J0[i], dpl[i], qpt[it][i],
								rid[it][i], hdcc0[it][i], vpt[it][i], Hwup[it][i], Hwdw[it][i], 
								ZJup[i], ZJdw[i],slop[i], Hj[I0[i]]);
				printStream.println();
			}
			// outfile<<endl;
			printStream.println();
			// -------------- ��ʼ���������ڵ� ---------------
			//
			for (i = 0; i < NP; i++) {
				j = I0[i];
				if (Hwup[it][i] == Hj[j]) {
					overflow[it][j] = overflow[it][j]
							+ (qpt[it][i] - qqkp[it][i]) * dt * 60.0;
					Hw_over[it][j] = csf * overflow[it][j] / Aj[j] / 10000.0
							* 1000.0;
				}
				if (Hwup[it][i] < Hj[j] && overflow[it][j] > 0.0) {
					overflow[it][j] = overflow[it][j] * 0.90;
					Hw_over[it][j] = csf * overflow[it][j] / Aj[j] / 10000.0
							* 1000.0;
				}
			}
			// ------------------ ���������ڵ���� ---------------
			//
		}// 1-- it end ---
			//
			// ----------------��Ļ����������------
			// cout<<"------ ģ�ͼ���ȫ����� ------"<<endl;
		System.out.println("------ ģ�ͼ���ȫ����� ------");
		//
		// --------------------------------- ����ܶγ����ȼ����� ---------------
		// outfile<<" ======== ʱ�ιܶγ����� ========"<<endl;
		printStream.println(" ======== ʱ�ιܶγ����� ========");
		// outfile<<"  i=    ";
		printStream.println("  i=    ");
		for (i = 0; i < NP; i++) {
			if (i < 10) {// outfile<<"    "<<i<<"   ";
				printStream.print("    " + i + "   ");
			} else {// outfile<<"   "<<i<<"   ";
				printStream.print("   " + i + "   ");
			}
		} // outfile<<endl;
		printStream.println();
		// outfile<<"it="<<endl;
		printStream.println("it=");
		for (it = 0; it < NT; it++) {
			if (it < 10) {// outfile<<" "<<it<<"   ";
				printStream.print(" " + it + "   ");
			} else {// outfile<<it<<"   ";
				printStream.print(it + "   ");
			}
			//
			for (i = 0; i < NP; i++) { // s.Format("%8.3f",hdcc0[it][i]);
				// outfile<<s;
				printStream.printf("%8.3f", hdcc0[it][i]);
			}
			// outfile<<endl;
			printStream.println();
		}
		//
		// --------------------------------- ����ڵ����������� ---------------
		// outfile<<" ======== ʱ�νڵ��ˮ��(m3) ========"<<endl;
		printStream.println(" ======== ʱ�νڵ��ˮ��(m3) ========");
		// outfile<<"  i=    ";
		printStream.print("  i=    ");
		for (i = 0; i < NN; i++) {
			if (i < 10) {// outfile<<"    "<<i<<"     ";
				printStream.print("    " + i + "     ");
			} else {// outfile<<"   "<<i<<"     ";
				printStream.print("   " + i + "     ");
			}
		} // outfile<<endl;
		printStream.println();
		// outfile<<"it="<<endl;
		printStream.println("it=");
		for (it = 0; it < NT; it++) {
			if (it < 10) {// outfile<<" "<<it<<"   ";
				printStream.print(" " + it + "   ");
			} else {// outfile<<it<<"   ";
				printStream.print(it + "   ");
			}
			//
			for (i = 0; i < NN; i++) {
				if (overflow[it][i] <= 0.0) {// outfile<<"          ";
					printStream.print("          ");
				} else { // s.Format("%10.2f",overflow[it][i]);
							// outfile<<s;
					printStream.printf("%10.2f", overflow[it][i]);
				}
			}
			// outfile<<endl;
			printStream.println();
		}
		//
		// outfile<<" ======== ʱ�νڵ��ˮ���(mm) ========"<<endl;
		printStream.println(" ======== ʱ�νڵ��ˮ���(mm) ========");
		// outfile<<"  i=    ";
		printStream.print("  i=    ");
		for (i = 0; i < NN; i++) {
			if (i < 10) {// outfile<<"    "<<i<<"     ";
				printStream.print("    " + i + "     ");
			} else {// outfile<<"   "<<i<<"     ";
				printStream.print("   " + i + "     ");
			}
		} // outfile<<endl;
		printStream.println();
		// outfile<<"it="<<endl;
		printStream.println("it=");
		for (it = 0; it < NT; it++) {
			if (it < 10) {// outfile<<" "<<it<<"   ";
				printStream.print(" " + it + "   ");
			} else {// outfile<<it<<"   ";
				printStream.print(it + "   ");
			}
			//
			for (i = 0; i < NN; i++) {
				if (overflow[it][i] <= 0.0) {// outfile<<"          ";
					printStream.print("          ");
				} else { // s.Format("%10.2f",Hw_over[it][i]);
							// outfile<<s;
					printStream.printf("%10.2f", Hw_over[it][i]);
				}
			}
			// outfile<<endl;
			printStream.println();
		}
	}
	// ====================================================================
	// outfile.close();
}
