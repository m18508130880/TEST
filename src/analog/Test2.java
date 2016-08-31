package analog;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.text.DecimalFormat;

public class Test2 {

	// 20150808new��ˮ��������ģ��Ƽ������-֥�Ӹ������
	// ����SQ LIU, TONGJI UNIVERSITY, 8 AUG 2016
	public static void main(String[] args) throws FileNotFoundException {

		String FileName = "201608031-��-��ˮ��������ģ��-֥�Ӹ������-���ҳ�002.txt";
		FileOutputStream fs = new FileOutputStream(new File(FileName));
		PrintStream printStream = new PrintStream(fs);
		printStream.println(FileName);

		DecimalFormat df = new DecimalFormat("##.######");

		// �����������ݣ���ϵͳ�������ܶ��������ڵ��������ܵ��������·�����ܶ�����������������м���д�ļ�ָ��
		// int Nbr = 8;//��ϵͳ��
		int NP = 9;// �ܶ���
		int NN = 10;// �ڵ���
		int Nstart = 3;// �ܵ������(�Զ�)
		int Npline = 7;// ·�����ܶ���(�Զ�)
		int Nmax = 20;// ���������
		int Iprt = 0;// �м���д�ļ�ָ��
		// ģ��ʱ������֥�Ӹ���ʱ��λ�ã�r=0.375�����ܵ�·������·�����ڵ������յ�ڵ��
		int NT = 60;// ģ��ʱ����
		int NR = 23;// ֥�Ӹ���ʱ��λ�ã�r=0.375��
		int Nroute = 3;// �ܵ�·���� (�Զ�)
		int Nr_node = 8;// ·�����ڵ���(�Զ�)
		int Nend = 8;// �յ�ڵ��

		/**
		 * �ܾ����ݼ��� ��ˮ���Aj ����ϵ��Acoef ������ Hj ���ױ��Hj ������� �ײ���� ֱ�� YJ002001 5.244
		 * 3.894 0.7 YG002000, YG002001 �� �� ���� ԭʼ̽�� ���ҳ�У����ʾ YJ002002 5.191 3.842
		 * 0.7 YG002001, YG002002 �м�� ���� ԭʼ̽�� ���ҳ�У����ʾ YJ002003 5.177 3.784 0.7
		 * YG002002, YG002003 �м�� ���� ԭʼ̽�� ���ҳ�У����ʾ YJ002004 5.208 3.733 0.7
		 * YG002003, YG002004 �м�� ���� ԭʼ̽�� ���ҳ�У����ʾ YJ002005 5.221 3.687 0.7
		 * YG002004, YG002005 �м�� ���� ԭʼ̽�� ���ҳ�У����ʾ YJ002006 5.201 3.643 0.7
		 * YG002005, YG002006 �м�� ���� ԭʼ̽�� ���ҳ�У����ʾ YJ002007 5.2 3.25 0.7
		 * YG002006,YG002009 YG002007 �м�� ���� ԭʼ̽�� ���ҳ�У����ʾ YJ002008 5.121 3.171
		 * 0.7 YG002007,YG002008 YG002999 �� �� ���� ԭʼ̽�� ���ҳ�У����ʾ YJ002009 5.131
		 * 3.731 0.7 YG002000, YG002008 �� �� ���� ԭʼ̽�� ���ҳ�У����ʾ YJ002010 5.186 3.486
		 * 0.7 YG002000, YG002009 �� �� ���� ԭʼ̽�� ���ҳ�У����ʾ
		 */
		// ������ʼ�ڵ�ź���ʼ�ڵ�ܵ�����<m>
		// int[] NJstart = new int[]{0, 8, 9};
		// double[] HJstart = new double[]{3.894, 3.731, 3.486};
		// �ڵ��ˮ���(ha)3.5
		double Aj[] = new double[] { 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5,
				3.5, 3.5 };
		// �ڵ��ˮ�������ϵ��0.6
		double[] Acoef = new double[] { 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6,
				0.6, 0.6 };
		// �ڵ�����ߣ�m��[NN=23]
		double[] Hj = new double[] { 5.244, 5.191, 5.177, 5.208, 5.221, 5.201,
				5.2, 5.121, 5.131, 5.186 };

		/**
		 * �ܶ����ݼ��� ���ܾ�I0 �յ�ܾ�J0 ����lp ֱ��dpl Ħ��ϵ��slp ��˱��ZJup �ն˱��ZJdw ֱ�� ���� ��˱��
		 * �ն˱�� YG002001 0.3 28.5 YJ002001 YJ002002 3.894 3.842 PE ԭʼ̽�� ���ҳ�У����ʾ
		 * YG002002 0.3 32 YJ002002 YJ002003 3.842 3.784 PE ԭʼ̽�� ���ҳ�У����ʾ
		 * YG002003 0.3 28.6 YJ002003 YJ002004 3.784 3.733 PE ԭʼ̽�� ���ҳ�У����ʾ
		 * YG002004 0.3 25.4 YJ002004 YJ002005 3.733 3.687 PE ԭʼ̽�� ���ҳ�У����ʾ
		 * YG002005 0.3 24.7 YJ002005 YJ002006 3.687 3.643 PE ԭʼ̽�� ���ҳ�У����ʾ
		 * YG002006 0.3 23.5 YJ002006 YJ002007 3.643 3.601 PE ԭʼ̽�� ���ҳ�У����ʾ
		 * YG002007 0.3 30.4 YJ002007 YJ002008 3.601 3.546 PE ԭʼ̽�� ���ҳ�У����ʾ
		 * YG002008 0.3 15.5 YJ002009 YJ002008 3.731 3.171 PE ԭʼ̽�� ���ҳ�У����ʾ
		 * YG002009 0.3 4.3 YJ002010 YJ002007 3.886 3.7 PE ԭʼ̽�� ���ҳ�У����ʾ
		 */
		// �ܶ����νڵ��I0,���νڵ��J0,
		int[] I0 = new int[] { 0, 1, 2, 3, 4, 5, 6, 8, 9 };
		int[] J0 = new int[] { 1, 2, 3, 4, 5, 6, 7, 7, 6 };
		// �ܶγ���
		double[] lp = new double[] { 28.5, 32, 28.6, 25.4, 24.7, 23.5, 30.4,
				15.5, 4.3 };
		// �ܶ�ֱ��(m)
		double[] dpl = new double[] { 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3,
				0.3 };
		// Ħ|��ϵ��
		double[] slp = new double[] { 0.017, 0.017, 0.017, 0.017, 0.017, 0.017,
				0.017, 0.017, 0.017 };
		// ���ιܵ׸߳�(m)�����ιܵ׸߳�(m)
		double[] ZJup = new double[] { 3.894, 3.842, 3.784, 3.733, 3.687,
				3.643, 3.601, 3.731, 3.886 };
		double[] ZJdw = new double[] { 3.842, 3.784, 3.733, 3.687, 3.643,
				3.601, 3.546, 3.171, 3.7 };
		/**
		 * �ڵ���� 1 2 3 4 5 6 7 8 
		 * ·��1	 0 1 2 3 4 5 6 7 
		 * ·��2	 8 7 -99 -99 -99 -99 -99 -99
		 * ·��3 	 9 6 -99 -99 -99 -99 -99 -99
		 */
		// ����·������·���ڵ��(-99��ʾ�սڵ�)
		int[][] Mroute = new int[][] { { 0, 1, 2, 3, 4, 5, 6, 7 },
				{ 8, 7, -99, -99, -99, -99, -99, -99 }, { 9, 6, -99, -99, -99, -99, -99, -99 } };

		/**
		 * �ڵ���� 1 2 3 4 5 6 7 ·��1 7 6 5 4 3 2 1 ·��2 8 -99 ·��3 9 -99
		 */
		// ��ϵͳ��֧·���ܶ����ݾ��� ����pipe branches-reverse order
		int[][] Mbranch = new int[][] { { 6, 5, 4, 3, 2, 1, 0 }, { 7, -99, -99, -99, -99, -99, -99 },
				{ 8, -99, -99, -99, -99, -99, -99 } };

		// double
		// dtnt,taa,tbb,AA,XX1,XX2,XX[NT],qit[NT],sumqj[NT][NN],sumAj[NT][NN],Tnode[NN][NN],sumTnode[NN][NN];
		double dtnt, taa, tbb, AA, XX1, XX2;
		double[] XX = new double[NT];
		double[] qit = new double[NT];
		double[][] sumqj = new double[NT][NN];
		double[][] sumAj = new double[NT][NN];
		double[][] Tnode = new double[NN][NN];
		double[][] sumTnode = new double[NN][NN];
		int i, j, k, ik, jk, k00, kkk1;
		int ip, it, k1, kp, in1, in2, in3;

		double hdc_min, hdc_max;
		double[] hdcc = new double[2];
		double[] fun_hd = new double[2];
		double fun_hd0 = 0;
		double[][] hdcc0 = new double[NT][NP];
		double hda, hdb, hdf, hdd, hdc = 0, hde = 0;
		// �ܶ����٣�m/s��
		double[] vp = new double[NP];
		double[] slop = new double[NP];
		double Ad0, AD, RD, vp0 = 1.0;
		double[][] qpt = new double[NT][NP];
		double[][] vpt = new double[NT][NP];
		double[][] rid = new double[NT][NP];
		double[][] slopt = new double[NT][NP];
		double[][] Hwup = new double[NT][NP];
		double[][] Hwdw = new double[NT][NP];
		double[][] qjt = new double[NT][NN];
		double[][] overflow = new double[NT][NN];
		double[][] Hw_over = new double[NT][NN];
		double hdj0, ARD23;
		double[] totalflow = new double[NN];
		double[] totalHw = new double[NN];
		// ���깫ʽshanghai storm water formular:
		// (A1+C*lgP)/(t+b)**n---ln(N)=2.303log(N)---����ˮλ��m��
		double A1 = 17.53, C_storm = 0.95, tmin = 10, b_storm = 11.77, P_simu = 10, n_storm = 0.88, dt = 2.0, rc = 0.375, Hw_end = 58.9;

		// --��������ļ���ʼ---
		// ofstream outfile;
		// outfile.open("20160808-��ˮ��������ģ��-֥�Ӹ������10-60.txt");
		//
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
		for (i = 0; i < NP; i++) { // s.Format("%4d",i);
									// outfile<<s;
			printStream.printf("%4d", i);
		}
		// outfile<<endl;
		printStream.println();

		// outfile<<" I0=";
		printStream.print(" I0=");
		for (i = 0; i < NP; i++) {// s.Format("%4d",I0[i]);
									// outfile<<s;
			printStream.printf("%4d", I0[i]);
		}
		// outfile<<endl;
		printStream.println();
		// outfile<<" J0=";
		printStream.print(" J0=");
		for (i = 0; i < NP; i++) {// s.Format("%4d",J0[i]);
									// outfile<<s;
			printStream.printf("%4d", J0[i]);
		}
		// outfile<<endl;
		// outfile<<endl;
		printStream.println();
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
		printStream.print("      ");
		for (j = 0; j < NN; j++) { // s.Format("%6d",j);
									// outfile<<s;
			printStream.printf("%6d", j);
		}
		// outfile<<endl;
		printStream.println();
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
		printStream.print("==j=  ");
		// outfile<<"      ";
		for (j = 0; j < NN; j++) { // s.Format("%6d",j);
									// outfile<<s;
			printStream.printf("%6d", j);
		}
		// outfile<<endl;
		printStream.println();
		for (i = 0; i < NN; i++) { // outfile<<"i="<<i<<"   ";
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
		// ================= ������̬����ģ��============================
		//
		// -----------��̬ģ����������------------------------------------
		// ----------------�ڵ��ˮ���(ha)�ͻ�ˮ����(m3/sec)����--------
		// outfile<<endl;
		printStream.println();
		// outfile<<"===========  ������̬ģ�����      �����ڣ� "<<P_simu<<"  ��   ʱ������ "<<NT<<"       �յ�ˮλ�� "<<Hw_end<<"  m  ========="<<endl;
		printStream.println("===========  ������̬ģ�����      �����ڣ� " + P_simu
				+ "  ��   ʱ������ " + NT + "       �յ�ˮλ�� " + Hw_end
				+ "  m  =========");
		// ----------rainfall intensity at every time step-----------
		// xxxxxxx
		// ֥�Ӹ������
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
			for (j = 0; j < NP; j++)
				qpt[it][j] = -99.0;
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
			// -------------------????????????????????????????????090127?????????????????------------------------
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
							printStream.println("   it= " + it + "   kp= " + kp
									+ "  Hwdm= " + df.format(Hwdw[it][kp])
									+ "  Hw_end= " + Hw_end);
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
							slopt[it][kp] = 10.29 * Math.pow(slp[kp], 2.0)
									* Math.pow(qpt[it][kp], 2.0)
									/ Math.pow(dpl[kp], 5.333);
							Hwup[it][kp] = Hwdw[it][kp] + slopt[it][kp]
									* lp[kp];
							if (Hwup[it][kp] >= Hj[I0[kp]]) {
								Hwup[it][kp] = Hj[I0[kp]];
								slopt[it][kp] = (Hwup[it][kp] - Hwdw[it][kp])
										/ lp[kp];
								vpt[it][kp] = Math.pow(rid[it][kp], 0.6667)
										* Math.pow(slopt[it][kp], 0.5)
										/ slp[kp];
								qpt[it][kp] = vpt[it][kp] * Ad0;
							}
						} else
						// --5
						{
							hdc_min = (Hwdw[it][kp] - ZJdw[kp]) / dpl[kp];
							slopt[it][kp] = slop[kp];
							if (hdc_min < 0.0)
								hdc_min = 0.0;
							ARD23 = slp[kp] * qpt[it][kp]
									/ Math.pow(slopt[it][kp], 0.5);
							// ----------h/d---0.618aaa---------------------------
							hdc_max = 1.0;
							fun_hd0 = 1000.0;
							for (kkk1 = 1; kkk1 < Nmax && fun_hd0 > 0.0001; kkk1++)
							// --6
							{
								if (kkk1 == 1) {
									hdcc[0] = hdc_min;
									hdcc[1] = 1.0;
								} else {
									hdcc[0] = hdc_min + (hdc_max - hdc_min)
											* 0.382;
									hdcc[1] = hdc_min + (hdc_max - hdc_min)
											* 0.618;
								}
								for (k00 = 0; k00 < 2; k00++) {
									hda = 1.0 - 2 * hdcc[k00];
									hdb = 1.0 - hdcc[k00];
									hdc = Math.pow(Math.cos(hda), -1);
									hdd = hdcc[k00] * hdb;
									hde = hda * Math.pow(hdd, 0.5);
									hdf = Math.pow(dpl[kp], 2.0);
									AD = hdf / 4.0 * hdc - hdf / 2.00 * hde;
									RD = dpl[kp] / 4.0 - dpl[kp] * hde
											/ (2.0 * hdc);
									fun_hd[k00] = AD * Math.pow(RD, 2.0 / 3.0)
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
								if (Iprt == 1)// outfile<<"it=  "<<it<<"  ik=  "<<ik<<"jk=  "<<jk<<"  kp=  "<<kp<<"  kkk1=  "<<kkk1<<"  hdcc[0]=  "<<hdcc[0]<<"  hdcc[1]=  "<<hdcc[1]
												// <<"  fun_hd[0]= "<<fun_hd[0]<<"  fun_hd[1]= "<<fun_hd[1]<<"  ARD23= "<<ARD23<<endl;
									printStream.println("it=  " + it
											+ "  ik=  " + ik + "jk=  " + jk
											+ "  kp=  " + kp + "  kkk1=  "
											+ kkk1 + "  hdcc[0]=  " + hdcc[0]
											+ "  hdcc[1]=  " + hdcc[1]
											+ "  fun_hd[0]= " + fun_hd[0]
											+ "  fun_hd[1]= " + fun_hd[1]
											+ "  ARD23= " + ARD23);
							} // 6-kkk1
								//
							rid[it][kp] = (0.25 - 0.5 * hde / hdc) * dpl[kp];
							Hwdw[it][kp] = ZJdw[kp] + hdcc0[it][kp] * dpl[kp];
							vpt[it][kp] = Math.pow(rid[it][kp], 2.0 / 3.0)
									* Math.pow(slopt[it][kp], 0.5) / slp[kp];
							Hwup[it][kp] = Hwdw[it][kp] + slopt[it][kp]
									* lp[kp];
							// ---------h/d---0.618aaa---------------------------
						} // 5--end_hdcc0[it][kp]
							//
							// outfile<<"   it= "<<it<<"   kp= "<<kp<<"   I0[kp]= "<<I0[kp]<<"  Hwdm= "<<Hwdw[it][kp]<<"  Hj= "<<Hj[I0[kp]]
							// <<"  hdcc0= "<<hdcc0[it][kp]<<"  qpt= "<<qpt[it][kp]<<"  vpt= "<<vpt[it][kp]<<"  hdcc0= "<<hdcc0[it][kp]<<" fun_hd0= "<<fun_hd0<<endl;
						printStream.println("   it= " + it + "   kp= " + kp
								+ "   I0[kp]= " + I0[kp] + "  Hwdm= "
								+ df.format(Hwdw[it][kp]) + "  Hj= "
								+ Hj[I0[kp]] + "  hdcc0= "
								+ df.format(hdcc0[it][kp]) + "  qpt= "
								+ df.format(qpt[it][kp]) + "  vpt= "
								+ df.format(vpt[it][kp]) + "  hdcc0= "
								+ df.format(hdcc0[it][kp]) + " fun_hd0= "
								+ df.format(fun_hd0));
					}// 4 if(kp>=0) end
				}// 3 ---jk end ---
			}// --2---ik end ---Hj[kp]
				// --------------------------------- ��ʼ���������ڵ� ---------------
				//
			for (i = 0; i < NP; i++) {
				j = I0[i];
				if (Hwup[it][i] == Hj[j]) {
					overflow[it][j] = overflow[it][j] - qpt[it][i];
					for (ip = 0; ip < NP; ip++) {
						k1 = J0[ip];
						if (k1 == j)
							overflow[it][j] = overflow[it][j] + qpt[it][ip];
					}
					qjt[it][j] = Aj[j] * qit[it] * Acoef[j];
					overflow[it][j] = (overflow[it][j] + qjt[it][j]) * dt
							* 60.0;
					Hw_over[it][j] = overflow[it][j] / Aj[j] / 10000.0 * 1000.0;
				}
			}
			// ----qjt[it][j],overflow[it][j],totalflow[j],totalHw[j]
			// ------------------ ���������ڵ���� ---------------
			// outfile<<endl;
			printStream.println();
			// outfile<<"    it   �ܶκ�  I0   J0   �ܾ�dpl   �ܶ�qp ˮ���뾶R    ������ ����(m/s)  ����ˮλ  ����ˮλ  �Ϲܵ׸�  �¹ܵ׸�  �ϵ����"<<endl;
			printStream
					.println("    it   �ܶκ�  I0   J0   �ܾ�dpl   �ܶ�qp ˮ���뾶R    ������ ����(m/s)  ����ˮλ  ����ˮλ  �Ϲܵ׸�  �¹ܵ׸�  �ϵ����");
			for (i = 0; i < NP; i++) {// s.Format("%6d%6d%6d%5d%10.2lf%10.3lf%10.3lf%10.3lf%10.3lf%10.3lf%10.3lf%10.3lf%10.3lf%10.3lf\n",
										// it,i,I0[i],J0[i],dpl[i],qpt[it][i],rid[it][i],hdcc0[it][i],vpt[it][i],Hwup[it][i],Hwdw[it][i],ZJup[i],ZJdw[i],Hj[I0[i]]);
				// outfile<<s;
				printStream
						.printf("%6d%6d%6d%5d%10.2f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f",
								it, i, I0[i], J0[i], dpl[i], qpt[it][i],
								rid[it][i], hdcc0[it][i], vpt[it][i],
								Hwup[it][i], Hwdw[it][i], ZJup[i], ZJdw[i],
								Hj[I0[i]]);
				printStream.println();
			}
			// outfile<<endl;
			printStream.println();
		}// 1-- it end ---
			// ------- ����ڵ��ˮ���ͻ�ˮ���(m)----
		for (j = 0; j < NN; j++) {
			totalflow[j] = 0.0;
			totalHw[j] = 0.0;
		}
		for (j = 0; j < NN; j++) {
			for (it = 0; it < NT; it++) {
				totalflow[j] = totalflow[j] + overflow[it][j];
				totalHw[j] = totalHw[j] + Hw_over[it][j];
			} // outfile<<"  totalflow[j]="<<totalflow[j]<<"   totalHw[j]="<<totalHw[j]<<endl;
		}
		// -----��Ļ����ܶ�ˮ���������------
		// cout<<"    it   �ܶκ�  I0   J0   �ܾ�dpl   �ܶ�qp ˮ���뾶R    ������ ����(m/s)  ����ˮλ  ����ˮλ  �Ϲܵ׸�  �¹ܵ׸�  �ϵ����"<<endl;
		System.out
				.println("    it   �ܶκ�  I0   J0   �ܾ�dpl   �ܶ�qp ˮ���뾶R    ������ ����(m/s)  ����ˮλ  ����ˮλ  �Ϲܵ׸�  �¹ܵ׸�  �ϵ����");
		// outfile<<endl;
		for (it = 0; it < NT; it++) {// outfile<<"  it= "<<it<<endl;
										// outfile<<endl;
			// outfile<<"    it   �ܶκ�  I0   J0   �ܾ�dpl   �ܶ�qp ˮ���뾶R    ������ ����(m/s)  ����ˮλ  ����ˮλ    �Ϲܵ׸߳� �¹ܵ׸߳�"<<endl;
			for (i = 0; i < NP; i++) {// s.Format("%6d%6d%6d%5d%10.2lf%10.3lf%10.3lf%10.3lf%10.3lf%10.3lf%10.3lf%10.3lf%10.3lf%10.3lf\n",
										// it,i,I0[i],J0[i],dpl[i],qpt[it][i],rid[it][i],hdcc0[it][i],vpt[it][i],Hwup[it][i],Hwdw[it][i],ZJup[i],ZJdw[i],Hj[I0[i]]);
										// outfile<<s;
				// cout<<s;
				System.out
						.printf("%6d%6d%6d%5d%10.2f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f",
								it, i, I0[i], J0[i], dpl[i], qpt[it][i],
								rid[it][i], hdcc0[it][i], vpt[it][i],
								Hwup[it][i], Hwdw[it][i], ZJup[i], ZJdw[i],
								Hj[I0[i]]);
				System.out.println();
			}
			// outfile<<endl;
		}
		// outfile<<endl;
		// --------------------------------- �ڵ����������� ---------------
		// outfile<<" ======== ʱ�νڵ��ˮ��(m3) ========"<<endl;
		printStream.println(" ======== ʱ�νڵ��ˮ��(m3) ========");

		// outfile<<"  i=    ";
		printStream.print("  i=    ");
		for (i = 0; i < NN; i++) {
			if (i < 10) {// outfile<<"  "<<i<<"   ";
				printStream.print("  " + i + "   ");
			} else {// outfile<<" "<<i<<"   ";
				printStream.print(" " + i + "   ");
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
				if (overflow[it][i] <= 0.0) {// outfile<<"      ";
					printStream.print("   ");
				} else { // s.Format("%6.1f",overflow[it][i]);
							// outfile<<s;
					printStream.printf("%6.1f", overflow[it][i]);
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
			if (i < 10) {// outfile<<"  "<<i<<"   ";
				printStream.print("  " + i + "   ");
			} else {// outfile<<" "<<i<<"   ";
				printStream.print(" " + i + "   ");
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
				if (overflow[it][i] <= 0.0) {// outfile<<"      ";
					printStream.print("      ");
				} else {// s.Format("%6.1f",Hw_over[it][i]);
						// outfile<<s;
					printStream.printf("%6.1f", Hw_over[it][i]);
				}
			}
			// outfile<<endl;
			printStream.println();
		}
		// ====================================================================
		// outfile.close();
	}

}
