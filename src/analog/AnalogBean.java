package analog;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;
import java.text.DecimalFormat;
import java.util.ArrayList;

import javax.servlet.ServletConfig;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;

import jxl.Sheet;
import jxl.Workbook;
import jxl.read.biff.BiffException;
import rmi.Rmi;
import util.CommUtil;
import util.CurrStatus;

import com.jspsmart.upload.SmartUpload;

public class AnalogBean
{
	
	/**
	 * 模拟计算时导入excel表格
	 * @param request
	 * @param response
	 * @param pRmi
	 * @param pFromZone
	 * @param pConfig
	 */
	public void ImportData(HttpServletRequest request, HttpServletResponse response, Rmi pRmi, boolean pFromZone, ServletConfig pConfig)
	{
		SmartUpload mySmartUpload = new SmartUpload();
		try
		{
			mySmartUpload.initialize(pConfig, request, response);
			mySmartUpload.setAllowedFilesList("xls,xlsx,XLS,XLSX,");
			mySmartUpload.upload();

			this.Sid = mySmartUpload.getRequest().getParameter("Sid");
			CurrStatus currStatus = (CurrStatus) request.getSession().getAttribute("CurrStatus_" + this.Sid);
			currStatus.getHtmlData(request, pFromZone);
			String Project_Id = mySmartUpload.getRequest().getParameter("Project_Id");

			if ((mySmartUpload.getFiles().getCount() > 0) && (mySmartUpload.getFiles().getFile(0).getFilePathName().trim().length() > 0))
			{
				if (mySmartUpload.getFiles().getFile(0).getSize() / 1024 <= 3072)
				{
					this.FileSaveRoute = "/www/DPP-LOCAL/DPP-LOCAL-WEB/files/analogData/";

					com.jspsmart.upload.File myFile = mySmartUpload.getFiles().getFile(0);
					this.File_Name = mySmartUpload.getFiles().getFile(0).getFileName();
					myFile.saveAs(this.FileSaveRoute + Project_Id + "_" + this.File_Name);
				}
				else
				{
					currStatus.setResult("文档上传失败！文档过大，必须小于3M!");
				}
			}
			currStatus.setJsp("loading.jsp?Sid=" + this.Sid + "&Project_Id=" + Project_Id);
			request.getSession().setAttribute("CurrStatus_" + this.Sid, currStatus);
			response.sendRedirect(currStatus.getJsp());
		}
		catch (Exception e)
		{
			e.printStackTrace();
		}
	}

	/**
	 * 计算管井时段水位 - 水位折綫圖
	 * @param gjId
	 * @return WaterAccGj
	 */
	public String AnalogWaterAccGj(String gjId)// subSys, timePeriod
	{
		AnalogWaterType = "WaterAccGj";
		return analog2(null, 0, gjId, AnalogWaterType);
	}

	/**
	 * 计算时段水位深度 - 水位剖面圖
	 * @param subSys
	 * @param timePeriod
	 * @return WaterLev
	 */
	public String AnalogWaterLev(String subSys, int timePeriod)
	{
		AnalogWaterType = "WaterLev";
		return analog2(subSys, timePeriod, null, AnalogWaterType);
	}

	/**
	 * 计算时段积水量 - 模擬地圖點位積水量
	 * @param fileName
	 * @param timePeriod
	 * @return WaterAcc
	 */
	public String AnalogWaterAcc(String subSys)
	{
		AnalogWaterType = "WaterAcc";
		return analog2(subSys, 0, null, AnalogWaterType);
	}

	// 第一套版本
	private String analog1(String subSys, int timePeriod, String gjId, String AnalogWaterType)
	{
		int SubgjId = 0;
		if (gjId != null)
		{
			SubgjId = CommUtil.StrToInt(gjId.substring(12, 15)) - 1;
		}
		try
		{
			// 管网基础数据：
			// 管段数，节点数，管道起点数，路径最大管段数，最大计算次数，模拟时段数，芝加哥峰点时段位置
			// 管道路径数，路径最大节点数，终点节点号，中间结果输出文件指针
			int NP = 9, NN = 10, Nstart = 3, Npline = 7, NT = 60, NR = 23, Nroute = 3, Nr_node = 8, Nend = 7, Iprt = 0;
			// 暴雨公式参数shanghai storm water formular:
			// (A1+C*lgP)/(t+b)**n---ln(N)=2.303log(N)---出口水位（m）
			// 管段流速（m/s）, 管段设定流速vp0，地面凹凸系数csf
			double A1 = 17.53, C_storm = 0.95, tmin = 10, b_storm = 11.77, P_simu = 50, n_storm = 0.88, dt = 2.0, rc = 0.375, Hw_end = 3.0, vp0 = 0.8, csf = 3.0;

			// 子系统管段数据：
			int[] I0; // 管段上游节点号I0,
			int[] J0; // 下游节点号J0,
			double[] lp; // 管段长度
			double[] dpl; // 管段直径(m)
			double[] slp; // 摩|阻系数
			double[] ZJup; // 上游管底高程(m)
			double[] ZJdw; // 下游管底高程(m)

			// 子系统节点数据
			// 管网起始节点号和起始节点管底埋深<m>
			double[] Aj; // 节点汇水面积(ha)3.5
			double[] Acoef; // 节点汇水面积径流系数0.6
			double[] Hj; // 节点地面标高（m）[NN=23]

			// 管网路径数和路径节点号(-99表示空节点)
			int[][] Mroute;

			// 子系统分支路径管段数据矩阵 倒序pipe branches-reverse order
			int[][] Mbranch;

			this.FileSaveRoute = "/www/DPP-LOCAL/DPP-LOCAL-WEB/files/analogData/";
			String XlsPath = "";
			if (gjId != null)
			{
				XlsPath = FileSaveRoute + gjId.substring(0, 12) + ".xls";
			}
			else
			{
				XlsPath = FileSaveRoute + subSys + ".xls";
			}
			InputStream is = new FileInputStream(XlsPath);
			Workbook rwb = Workbook.getWorkbook(is);
			Sheet rs = rwb.getSheet(0);
			int rsRows = rs.getRows();

			/*
			 * 基础数据表格子系统号 节点数NN 管段数NP 起点数NStart 路径管段数Npline 路径节点数Nr_node
			 * 终点出口号Nend 模拟时段NT 管段路径数NrouteYJ002 10 9 3 7 8 8 60 3
			 */
			int rowCnt = 2;
			String sysName = rs.getCell(0, rowCnt).getContents().trim();
			NN = Integer.parseInt(rs.getCell(1, rowCnt).getContents().trim());
			NP = Integer.parseInt(rs.getCell(2, rowCnt).getContents().trim());
			Nstart = Integer.parseInt(rs.getCell(3, rowCnt).getContents().trim());
			Npline = Integer.parseInt(rs.getCell(4, rowCnt).getContents().trim());
			Nr_node = Integer.parseInt(rs.getCell(5, rowCnt).getContents().trim());
			Nend = Integer.parseInt(rs.getCell(6, rowCnt).getContents().trim());
			NT = Integer.parseInt(rs.getCell(7, rowCnt).getContents().trim());
			Nroute = Integer.parseInt(rs.getCell(8, rowCnt).getContents().trim());
			rowCnt += 4;

			/*
			 * 子系统管段数据表格 Pipe.No 起点号I0 终点号J0 长度LP 直径DP 摩阻系数 起端标高 终端标高 1 0 1 28.5
			 * 0.3 0.017 3.894 3.842 2 1 2 32 0.3 0.017 3.842 3.784 3 2 3 28.6
			 * 0.3 0.017 3.784 3.733 4 3 4 25.4 0.3 0.017 3.733 3.687 5 4 5 24.7
			 * 0.3 0.017 3.687 3.643 6 5 6 23.5 0.3 0.017 3.643 3.601 7 6 7 30.4
			 * 0.3 0.017 3.601 3.546 8 8 7 15.5 0.3 0.017 3.731 3.171 9 9 6 4.3
			 * 0.3 0.017 3.886 3.7
			 */
			I0 = new int[NP];
			J0 = new int[NP];
			lp = new double[NP];
			dpl = new double[NP];
			slp = new double[NP];
			ZJup = new double[NP];
			ZJdw = new double[NP];
			for (int j = 0; j < NP; j++)
			{
				I0[j] = Integer.parseInt(rs.getCell(1, rowCnt + j).getContents().trim());
				J0[j] = Integer.parseInt(rs.getCell(2, rowCnt + j).getContents().trim());
				lp[j] = Double.parseDouble(rs.getCell(3, rowCnt + j).getContents().trim());
				dpl[j] = Double.parseDouble(rs.getCell(4, rowCnt + j).getContents().trim());
				slp[j] = Double.parseDouble(rs.getCell(5, rowCnt + j).getContents().trim());
				ZJup[j] = Double.parseDouble(rs.getCell(6, rowCnt + j).getContents().trim());
				ZJdw[j] = Double.parseDouble(rs.getCell(7, rowCnt + j).getContents().trim());
			}
			rowCnt += NP;
			rowCnt += 3;

			/*
			 * 子系统节点数据表格节点No 汇水面积ha 径流系数 地面标高 井底标高 1 3.5 0.6 5.244 暂未用到 2 3.5
			 * 0.6 5.191 3 3.5 0.6 5.177 4 3.5 0.6 5.208 5 3.5 0.6 5.221 6 3.5
			 * 0.6 5.201 7 3.5 0.6 5.2 8 3.5 0.6 5.121 9 3.5 0.6 5.131 10 3.5
			 * 0.6 5.186
			 */
			Aj = new double[NN];
			Acoef = new double[NN];
			Hj = new double[NN];
			for (int j = 0; j < NN; j++)
			{
				Aj[j] = Double.parseDouble(rs.getCell(1, rowCnt + j).getContents().trim());
				Acoef[j] = Double.parseDouble(rs.getCell(2, rowCnt + j).getContents().trim());
				Hj[j] = Double.parseDouble(rs.getCell(3, rowCnt + j).getContents().trim());
			}
			rowCnt += NN;
			rowCnt += 3;

			/*
			 * 管网路径数&路径节点号节点序号 1 2 3 4 5 6 7 8 1 0 1 2 3 4 5 6 7 2 8 7 -99 -99
			 * -99 -99 -99 -99 3 9 6 -99 -99 -99 -99 -99 -99
			 */
			Mroute = new int[Nstart][Nr_node];
			for (int j = 0; j < Nstart; j++)
			{
				for (int k = 0; k < Nr_node; k++)
				{
					Mroute[j][k] = Integer.parseInt(rs.getCell(k + 1, rowCnt + j).getContents().trim());
				}
			}
			rowCnt += Nstart;
			rowCnt += 3;

			/*
			 * 子系统分支路径管段数据矩阵 倒序pipe branches-reverse order 节点序号 1 2 3 4 5 6 7 1
			 * 6 5 4 3 2 1 0 2 7 -99 -99 -99 -99 -99 -99 3 8 -99 -99 -99 -99 -99
			 * -99
			 */
			Mbranch = new int[Nstart][Npline];
			for (int j = 0; j < Nstart; j++)
			{
				for (int k = 0; k < Npline; k++)
				{
					Mbranch[j][k] = Integer.parseInt(rs.getCell(k + 1, rowCnt + j).getContents().trim());
				}
			}
			// ----临界水深计算变量----
			double sita0 = 3.0, eps = 0.001, alfa = 0.5;
			double Ad0, qkpmax, Hwdwkp, yykp, sita, cons_b, sita_s = 0, sita_c, fsita, dfdsita, dfsita, ssita = 0, csita = 0, hyd_A, hafsita, shafsita = 0, chafsita, sita_p = 0;
			// 中间变量
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
			double[][] Hwj = new double[NT][NN];

			// ----------------------------------------------------------------------------------------------------------
			String FileName = "";
			if (gjId != null)
			{
				FileName = gjId.substring(0, 12) + ".txt";
			}
			else
			{
				FileName = subSys + ".txt";
			}
			String FilePath = "./www/DPP-LOCAL/DPP-LOCAL-WEB/files/analogValue/";
			FileOutputStream fs = new FileOutputStream(new File(FilePath + FileName));
			PrintStream printStream = new PrintStream(fs);
			printStream.println(FileName);

			DecimalFormat df = new DecimalFormat("##.####");
			DecimalFormat df1 = new DecimalFormat("######.##");
			// --输出数据文件开始---
			// ================= 赋初值 ===============================
			for (i = 0; i < NT; i++)
			{
				for (j = 0; j < NN; j++)
					sumAj[i][j] = 0;
			}
			for (i = 0; i < NT; i++)
			{
				for (j = 0; j < NN; j++)
					sumqj[i][j] = 0;
			}
			for (i = 0; i < NN; i++)
			{
				for (j = 0; j < NN; j++)
				{
					if (i == j)
					{
						Tnode[i][j] = 0;
					}
					else
					{
						Tnode[i][j] = -99;
					}
				}
			}
			for (i = 0; i < NN; i++)
			{
				for (j = 0; j < NN; j++)
					sumTnode[i][j] = 0;
			}
			// ==================Tnode-sumTnode=========================
			for (i = 0; i < NP; i++)
				vp[i] = vp0;
			for (kp = 0; kp < NP; kp++)
			{
				in1 = I0[kp];
				in2 = J0[kp];
				Tnode[in1][in2] = lp[kp] / vp[kp] / 60;
				slop[kp] = (ZJup[kp] - ZJdw[kp]) / lp[kp];
			}
			//
			for (i = 0; i < Nroute; i++)
			{
				for (j = 0; j < Nr_node; j++)
				{
					in1 = Mroute[i][j];
					if (in1 >= 0)
					{
						for (k = j + 1; k < Nr_node; k++)
						{
							in2 = Mroute[i][k - 1];
							in3 = Mroute[i][k];
							if (in3 >= 0)
							{
								sumTnode[in1][in3] = sumTnode[in1][in2] + Tnode[in2][in3];
							}
						}
					}
				}
			}
			// =====print Mroute[i][j], Tnode, sumTnode,Mbranch[i][j]====
			// System.out.println("pipe no.  I0    J0");
			for (i = 0; i < NP; i++)
			{
				// System.out.printf("%6d%6d%6d", i, I0[i], J0[i]);
				// System.out.println();
			}
			printStream.println();
			printStream.print(" ip=");
			for (i = 0; i < NP; i++)
			{
				printStream.printf("%4d", i);
			}
			printStream.println();
			printStream.print(" I0=");
			for (i = 0; i < NP; i++)
			{
				printStream.printf("%4d", I0[i]);
			}
			printStream.println();
			printStream.print(" J0=");
			for (i = 0; i < NP; i++)
			{
				printStream.printf("%4d", J0[i]);
			}
			printStream.println();
			printStream.println();
			printStream.println("===========  print Mroute[i][j]");
			for (i = 0; i < Nroute; i++)
			{
				for (j = 0; j < Nr_node; j++)
				{
					printStream.printf("%6d", Mroute[i][j]);
				}
				printStream.println();
			}
			printStream.println();
			printStream.println("===========  print Mbranch[i][j]");
			for (i = 0; i < Nstart; i++)
			{
				for (j = 0; j < Npline; j++)
				{
					printStream.printf("%6d", Mbranch[i][j]);
				}
				printStream.println();
			}
			printStream.println("===========  print Tnode[i][j]");
			printStream.println("====j=  ");
			printStream.println("      ");
			for (j = 0; j < NN; j++)
			{
				printStream.printf("%6d", j);
			}
			printStream.println();
			for (i = 0; i < NN; i++)
			{
				if (i < 10)
				{
					printStream.print("i=" + i + "   ");
				}
				else
				{
					printStream.print("i=" + i + "  ");
				}
				for (j = 0; j < NN; j++)
				{
					if (Tnode[i][j] < 0.0)
					{
						printStream.print("      ");
					}
					else
					{
						printStream.printf("%6.2f", Tnode[i][j]);
					}
				}
				printStream.println();
			}
			printStream.println();
			printStream.println("===========  print sumTnode[i][j]");
			printStream.print("==j=  ");
			for (j = 0; j < NN; j++)
			{
				printStream.printf("%6d", j);
			}
			printStream.println();
			for (i = 0; i < NN; i++)
			{
				printStream.print("i=" + i + "   ");
				for (j = 0; j < NN; j++)
				{
					if (sumTnode[i][j] <= 0.0)
					{
						printStream.print("      ");
					}
					else
					{
						printStream.printf("%6.2f", sumTnode[i][j]);
					}
				}
				printStream.println();
			}
			// ================= 管网准稳态流动模拟============================
			// -------------------动态模拟流量计算-----------------------------
			// ----------------节点汇水面积(ha)和汇水流量(m3/sec)计算--------
			printStream.println();
			printStream.println("===========  管网动态模拟计算      重现期＝ " + P_simu + "  年   时段数＝ " + NT + "       终点水位＝ " + Hw_end + "  m  =========");
			// 芝加哥过程线--rainfall intensity at every time step--
			AA = A1 + A1 * C_storm * Math.log(P_simu) / 2.303;
			for (it = 0; it < NT; it++)
			{
				if (it <= NR)
				{
					dtnt = dt * (float) (it);
					tbb = dt * (float) (NR) - dtnt;
					XX1 = AA * ((1.0 - n_storm) * tbb / rc + b_storm);
					XX2 = Math.pow((tbb / rc + b_storm), (n_storm + 1.0));
				}
				else
				{
					dtnt = dt * (float) (it);
					taa = dtnt - dt * (float) (NR);
					XX1 = AA * ((1.0 - n_storm) * taa / (1.0 - rc) + b_storm);
					XX2 = Math.pow((taa / (1.0 - rc) + b_storm), (n_storm + 1.0));
				}
				XX[it] = XX1 / XX2;
				qit[it] = 167.0 * XX[it] / 1000.0;
			}
			printStream.println();
			printStream.println("    it      dtnt      XX[it]     qit[it]");
			for (it = 0; it < NT; it++)
			{
				dtnt = dt * (float) (it);
				printStream.printf("%6d%10.2f%12.6f%12.6f", it, dtnt, XX[it], qit[it]);
				printStream.println();
			}
			printStream.println();
			for (it = 0; it < NT; it++)
			{
				dtnt = dt + dt * (float) (it);
				for (j = 0; j < NN; j++)
				{
					sumAj[it][j] = Aj[j];
					sumqj[it][j] = Aj[j] * qit[it] * Acoef[j];
					for (i = 0; i < NN; i++)
					{
						if (sumTnode[i][j] > 0 && sumTnode[i][j] < dtnt)
						{
							sumAj[it][j] = sumAj[it][j] + Aj[i];
							sumqj[it][j] = sumqj[it][j] + Aj[i] * qit[it] * Acoef[i];
						}
					}
				}
			}
			printStream.println("  sumAj[it][j]=");
			for (it = 0; it < NT; it++)
			{
				for (j = 0; j < NN; j++)
				{
					printStream.printf("%8.2f", sumAj[it][j]);
				}
				printStream.println();
			}
			printStream.println();
			printStream.println("  sumqj[it][j]=");
			for (it = 0; it < NT; it++)
			{
				for (j = 0; j < NN; j++)
				{
					printStream.printf("%8.2f", sumqj[it][j]);
				}
				printStream.println();
			}
			printStream.println();
			// ---------------------------------------------------------------
			for (it = 0; it < NT; it++)
			{
				for (i = 0; i < NN; i++)
				{
					overflow[it][i] = 0.0;
					Hw_over[it][i] = 0.0;
				}
			}
			for (it = 0; it < NT; it++)
			{
				for (j = 0; j < NP; j++)
				{
					qpt[it][j] = -99.0;
					qqkp[it][j] = 0.0;
				}
			}
			// ---------------------------------------------------------------
			for (it = 0; it < NT; it++)
			{
				printStream.print(" it=" + it + "  qpt[it][k]=");
				for (j = 0; j < NN; j++)
				{
					for (k = 0; k < NP; k++)
					{
						if (I0[k] == j)
						{
							qpt[it][k] = sumqj[it][j];
							printStream.printf("%8.2f", qpt[it][k]);
						}
					}
				}
				printStream.println();
				for (ik = 0; ik < Nstart; ik++)
				{
					for (jk = 0; jk < Npline; jk++)
					{
						kp = Mbranch[ik][jk];
						if (kp >= 0)
						{
							if (J0[kp] == Nend)
							{
								Hwdw[it][kp] = Hw_end;
								if (1 == Iprt)
								{
									printStream.println("   it= " + it + "   kp= " + kp + "  Hwdm= " + Hwdw[it][kp] + "  Hw_end= " + Hw_end);
								}
							}
							else
							{
								for (k1 = 0; k1 < NP; k1++)
								{
									if (I0[k1] == J0[kp]) Hwdw[it][kp] = Hwup[it][k1];
								}
							}
							Ad0 = 0.7854 * Math.pow(dpl[kp], 2.0);
							hdj0 = ZJdw[kp] + dpl[kp];
							if (Hwdw[it][kp] >= hdj0)
							{
								if (1 == Iprt)
								{
									printStream.println("   it= " + it + "   kp= " + kp + "  Hwdm= " + df.format(Hwdw[it][kp]) + "  淹没出流 ");
								}
								hdcc0[it][kp] = 1.0;
								rid[it][kp] = dpl[kp] / 4.0;
								vpt[it][kp] = qpt[it][kp] / Ad0;
								slopt[it][kp] = 10.29 * Math.pow(slp[kp], 2.0) * Math.pow(qpt[it][kp], 2.0) / Math.pow(dpl[kp], 5.333);
								Hwup[it][kp] = Hwdw[it][kp] + slopt[it][kp] * lp[kp];
								if (Hwup[it][kp] >= Hj[I0[kp]])
								{
									Hwup[it][kp] = Hj[I0[kp]];
									slopt[it][kp] = (Hwup[it][kp] - Hwdw[it][kp]) / lp[kp];
									if (slopt[it][kp] < 0.0)
									{
										slopt[it][kp] = Math.abs(slopt[it][kp]);
									}
									vpt[it][kp] = Math.pow(rid[it][kp], 0.6667) * Math.pow(slopt[it][kp], 0.5) / slp[kp];
									qqkp[it][kp] = vpt[it][kp] * Ad0;
									if (qqkp[it][kp] < 0.0)
									{
										qqkp[it][kp] = Math.abs(qqkp[it][kp]);
									}
								}
							}
							else
							{
								if (1 == Iprt)
								{
									printStream.println("   it= " + it + "   kp= " + kp + "  Hwdw= " + df.format(Hwdw[it][kp]) + "  非淹没出流 ");
								}
								// --20160907修改开始---采用临界水深算法-----------------------
								qkpmax = 2.46 * Math.pow(dpl[kp], 2.5);
								if (qpt[it][kp] > qkpmax * 0.95)
								{
									if (1 == Iprt)
									{
										printStream.println("   it= " + it + "   kp= " + kp + "  qkpmax= " + qkpmax + "  非淹没满管出流 ");
									}
									Hwdw[it][kp] = ZJdw[kp] + dpl[kp] * 1.1;
									hdcc0[it][kp] = 1.0;
									rid[it][kp] = dpl[kp] / 4.0;
									vpt[it][kp] = qpt[it][kp] / Ad0;
									slopt[it][kp] = 10.29 * Math.pow(slp[kp], 2.0) * Math.pow(qpt[it][kp], 2.0) / Math.pow(dpl[kp], 5.333);
									Hwup[it][kp] = Hwdw[it][kp] + slopt[it][kp] * lp[kp];
									if (Hwup[it][kp] >= Hj[I0[kp]])
									{
										Hwup[it][kp] = Hj[I0[kp]];
										slopt[it][kp] = (Hwup[it][kp] - Hwdw[it][kp]) / lp[kp];
										if (slopt[it][kp] < 0.0)
										{
											slopt[it][kp] = Math.abs(slopt[it][kp]);
										}
										vpt[it][kp] = Math.pow(rid[it][kp], 0.6667) * Math.pow(slopt[it][kp], 0.5) / slp[kp];
										qqkp[it][kp] = vpt[it][kp] * Ad0;
										if (qqkp[it][kp] < 0.0)
										{
											qqkp[it][kp] = Math.abs(qqkp[it][kp]);
										}
									}
								}
								else
								{
									if (1 == Iprt)
									{
										printStream.println("   it= " + it + "   kp= " + kp + "  Hwdm= " + df.format(Hwdw[it][kp]) + "  非淹没非满管出流 ");
									}
									i = 0;
									sita = sita0;
									cons_b = 0.276843 * Math.pow(dpl[kp], 2.5) / qpt[it][kp];
									if (Iprt == 1)
									{
										printStream.println("   k= " + k + "   qpt[it][kp]= " + qpt[it][kp] + "   cons_b= " + cons_b);
									}
									while (true)
									{
										ssita = Math.sin(sita);
										csita = Math.cos(sita);
										hafsita = sita / 2.0;
										shafsita = Math.sin(hafsita);
										chafsita = Math.cos(hafsita);
										sita_s = sita - Math.sin(sita);
										sita_c = 1 - Math.cos(sita);
										sita_p = Math.pow((1.0 - chafsita), -0.5);
										fsita = cons_b * sita_s - sita_p;
										dfsita = Math.abs(fsita);
										if (dfsita < eps)
										{
											hdcc0[it][kp] = (1 - Math.cos(sita / 2)) / 2;
											rid[it][kp] = 0.25 * dpl[kp] * (sita - Math.sin(sita)) / sita;
											vpt[it][kp] = Math.pow(rid[it][kp], 0.6667) * Math.pow(slop[kp], 0.5) / slp[kp];
											break;
										}
										else
										{
											dfdsita = cons_b * (1.0 - csita) + 0.25 * Math.pow(sita_p, -1.0) * shafsita;
											sita = sita - alfa * fsita / dfdsita;
											if (Iprt == 1)
											{
												printStream.println("   i= " + i + "   sita= " + sita + "   ssita= " + ssita + "   csita= " + csita + "   fsita= " + fsita + "   dfdsita= " + dfdsita);
											}
											i = i + 1;
										}
									}
								}
								Hwdwkp = ZJdw[kp] + hdcc0[it][kp] * dpl[kp];
								if (Hwdwkp >= Hwdw[it][kp])
								{
									Hwdw[it][kp] = Hwdwkp;
								}
								if (Hwdwkp < Hwdw[it][kp])
								{
									yykp = Hwdw[it][kp] - ZJdw[kp];
									if (yykp > dpl[kp])
									{
										yykp = dpl[kp];
									}
									sita = 2.0 * Math.acos(1.0 - 2.0 * yykp / dpl[kp]);
									hdcc0[it][kp] = yykp / dpl[kp];
									rid[it][kp] = 0.25 * dpl[kp] * (sita - Math.sin(sita)) / sita;
									vpt[it][kp] = Math.pow(rid[it][kp], 0.6667) * Math.pow(slop[kp], 0.5) / slp[kp];
								}
								Hwup[it][kp] = Hwdw[it][kp] + slop[kp] * lp[kp];
							}
							// ------- 输出it计算结果 ----------
							if (Iprt == 1)
							{
								printStream.println("   it= " + it + "   kp= " + kp + "   I0[kp]= " + I0[kp] + "  Hwdm= " + Hwdw[it][kp] + "  Hwup= " + Hwup[it][kp] + "  Hj= " + Hj[I0[kp]] + "  hdcc0= " + hdcc0[it][kp] + "  qpt= " + qpt[it][kp] + "  qqkp= " + qqkp[it][kp] + "  vpt= " + vpt[it][kp]);
							}
						}
					}
				}
				printStream.println();

				printStream.println("    it   管段号  I0   J0 管径dpl     管段qp   水力半径R  充满度 流速(m/s)  上游水位  下游水位  上管底高  下管底高  管段坡度  上地面高");
				for (i = 0; i < NP; i++)
				{
					printStream.printf("%6d%6d%6d%5d%8.2f%12.3f%10.3f%8.3f%10.3f%10.3f%10.3f%10.3f%10.3f%10.5f%10.3f", it, i, I0[i], J0[i], dpl[i], qpt[it][i], rid[it][i], hdcc0[it][i], vpt[it][i], Hwup[it][i], Hwdw[it][i], ZJup[i], ZJdw[i], slop[i], Hj[I0[i]]);
					printStream.println();
				}
				printStream.println();
				// -------------- 开始计算溢流节点 ---------------
				for (i = 0; i < NP; i++)
				{
					k = J0[i];
					if (k == Nend)
					{
						Hwj[it][k] = Hwdw[it][i];
					}
					{
						j = I0[i];
						Hwj[it][j] = Hwup[it][i];
						if (Hwup[it][i] == Hj[j])
						{
							overflow[it][j] = overflow[it][j] + (qpt[it][i] - qqkp[it][i]) * dt * 60.0;
							Hw_over[it][j] = csf * overflow[it][j] / Aj[j] / 10000.0 * 1000.0;

						}
						if (Hwup[it][i] < Hj[j] && overflow[it][j] > 0.0)
						{
							overflow[it][j] = overflow[it - 1][j] * 0.90;
							Hw_over[it][j] = csf * overflow[it][j] / Aj[j] / 10000.0 * 1000.0;
						}
					}
					if (it > NR && Hw_over[it][j] <= 5.0)
					{
						overflow[it][j] = 0.0;
						Hw_over[it][j] = 0.0;
					}
				}
				// ------------------ 计算溢流节点结束 ---------------
			}
			// ----------------屏幕输出计算结束------
			// System.out.println("------ 模型计算全部完成 ------");
			// ---------------- 输出管段充满度计算结果 ---------------
			printStream.println(" ======== 时段管段充满度 ========");
			printStream.print("  i=    ");
			for (i = 0; i < NP; i++)
			{
				if (i < 10)
				{
					printStream.print("    " + i + "   ");
				}
				else
				{
					printStream.print("   " + i + "   ");
				}
			}
			printStream.println();
			printStream.println("it=");
			for (it = 0; it < NT; it++)
			{
				if (it < 10)
				{
					printStream.print(" " + it + "   ");
				}
				else
				{
					printStream.print(it + "   ");
				}
				for (i = 0; i < NP; i++)
				{
					printStream.printf("%8.3f", hdcc0[it][i]);
				}
				printStream.println();
			}
			// --------------- 输出节点水位计算结果 ---------------
			printStream.println(" ======== 时段节点水位 ========");
			printStream.println("  i=    ");
			for (i = 0; i < NN; i++)
			{
				if (i < 10)
				{
					printStream.print("    " + i + "     ");
				}
				else
				{
					printStream.print("   " + i + "     ");
				}
			}
			printStream.println("it=");
			for (it = 0; it < NT; it++)
			{
				if (it < 10)
				{
					printStream.print(" " + it + "   ");
				}
				else
				{
					printStream.print(it + "   ");
				}
				String WaterLevNew = "";
				for (i = 0; i < NN; i++)
				{
					printStream.printf("%10.3f", Hwj[it][i]);
					if (gjId != null && i == SubgjId)
					{
						WaterAccGj += df1.format(Hwj[it][i]) + "|";
					}
					WaterLevNew += df1.format(Hwj[it][i]) + "|";
				}
				printStream.println();
				WaterLev[it] = WaterLevNew;
			}
			// ---------------- 输出节点溢流计算结果 ---------------
			printStream.println(" ======== 时段节点积水量(m3) ========");
			printStream.print("  i=    ");
			for (i = 0; i < NN; i++)
			{
				if (i < 10)
				{
					printStream.print("    " + i + "     ");
				}
				else
				{
					printStream.print("   " + i + "     ");
				}
			}
			printStream.println();
			printStream.println("it=");

			for (it = 0; it < NT; it++)
			{
				if (it < 10)
				{
					printStream.print(" " + it + "   ");
				}
				else
				{
					printStream.print(it + "   ");
				}
				String WaterAccNew = "";
				for (i = 0; i < NN; i++)
				{
					if (overflow[it][i] <= 0.0)
					{
						printStream.print("          ");
						WaterAccNew += 0 + "|";
					}
					else
					{
						printStream.printf("%10.2f", overflow[it][i]);
						WaterAccNew += df1.format(overflow[it][i]) + "|";
					}
				}
				WaterAcc[it] = WaterAccNew;
				printStream.println();
			}
			printStream.println(" ======== 时段节点积水深度(mm) ========");
			printStream.print("  i=    ");
			for (i = 0; i < NN; i++)
			{
				if (i < 10)
				{
					printStream.print("    " + i + "     ");
				}
				else
				{
					printStream.print("   " + i + "     ");
				}
			}
			printStream.println();
			printStream.println("it=");
			for (it = 0; it < NT; it++)
			{
				if (it < 10)
				{
					printStream.print(" " + it + "   ");
				}
				else
				{
					printStream.print(it + "   ");
				}
				for (i = 0; i < NN; i++)
				{
					if (overflow[it][i] <= 0.0)
					{
						printStream.print("          ");
					}
					else
					{
						printStream.printf("%10.2f", Hw_over[it][i]);
					}
				}
				printStream.println();
			}
		}
		catch (Exception e)
		{
			e.printStackTrace();
		}
		if (AnalogWaterType.equals("WaterAccGj"))
		{
			return WaterAccGj;
		}
		else if (AnalogWaterType.equals("WaterAcc"))
		{
			String WaterAccList = "";
			for (int i = 0; i < WaterAcc.length; i++)
			{
				WaterAccList += subSys.substring(7, 12) + WaterAcc[i] + ";";
			}
			return WaterAccList;
		}
		else if (AnalogWaterType.equals("WaterLev"))
		{
			String WaterLevList = "";
			for (int i = 0; i < WaterLev.length; i++)
			{
				WaterLevList += subSys.substring(7, 12) + WaterLev[i] + ";";
			}
			return WaterLevList;
			// return WaterLev[timePeriod];
		}
		return "";
	}

	// 第二套版本
	private String analog2(String subSys, int timePeriod, String gjId, String AnalogWaterType)
	{
		int SubgjId = 0;
		if (gjId != null)
		{
			SubgjId = CommUtil.StrToInt(gjId.substring(12, 15)) - 1;
		}
		try
		{
			// 管网基础数据：
			// 管段数，节点数，管道起点数，路径最大管段数，最大计算次数，模拟时段数，芝加哥峰点时段位置
			// 管道路径数，路径最大节点数，终点节点号，中间结果输出文件指针
			int NP = 9, NN = 10, Nstart = 3, Npline = 7, NT = 60, NR = 23, Nroute = 3, Nr_node = 8, Nend = 7, Iprt = 0, Nprtc = 20;
			// 暴雨公式参数shanghai storm water formular:
			// (A1+C*lgP)/(t+b)**n---ln(N)=2.303log(N)---出口水位（m）
			// 管段流速（m/s）, 管段设定流速vp0，地面凹凸系数csf
			double A1 = 17.53, C_storm = 0.95, tmin = 10, b_storm = 11.77, P_simu = 100, n_storm = 0.88, dt = 2.0, rc = 0.375, Hw_end = 3.0, vp0 = 0.8, csf = 3.0;

			// 子系统管段数据：
			int[] I0; // 管段上游节点号I0,
			int[] J0; // 下游节点号J0,
			double[] lp; // 管段长度
			double[] dpl; // 管段直径(m)
			double[] slp; // 摩|阻系数
			double[] ZJup; // 上游管底高程(m)
			double[] ZJdw; // 下游管底高程(m)

			// 子系统节点数据
			// 管网起始节点号和起始节点管底埋深<m>
			double[] Aj; // 节点汇水面积(ha)3.5
			double[] Acoef; // 节点汇水面积径流系数0.6
			double[] Hj; // 节点地面标高（m）[NN=23]

			// 管网路径数和路径节点号(-99表示空节点)
			int[][] Mroute;

			// 子系统分支路径管段数据矩阵 倒序pipe branches-reverse order
			int[][] Mbranch;

			this.FileSaveRoute = "/www/DPP-LOCAL/DPP-LOCAL-WEB/files/analogData/";
			String XlsPath = "";
			if (gjId != null)
			{
				XlsPath = FileSaveRoute + gjId.substring(0, 12) + ".xls";
			}
			else
			{
				XlsPath = FileSaveRoute + subSys + ".xls";
			}
			InputStream is = new FileInputStream(XlsPath);
			Workbook rwb = Workbook.getWorkbook(is);
			Sheet rs = rwb.getSheet(0);
			int rsRows = rs.getRows();

			/*
			 * 基础数据表格子系统号 节点数NN 管段数NP 起点数NStart 路径管段数Npline 路径节点数Nr_node
			 * 终点出口号Nend 模拟时段NT 管段路径数NrouteYJ002 10 9 3 7 8 8 60 3
			 */
			int rowCnt = 2;
			String sysName = rs.getCell(0, rowCnt).getContents().trim();
			NN = Integer.parseInt(rs.getCell(1, rowCnt).getContents().trim());
			NP = Integer.parseInt(rs.getCell(2, rowCnt).getContents().trim());
			Nstart = Integer.parseInt(rs.getCell(3, rowCnt).getContents().trim());
			Npline = Integer.parseInt(rs.getCell(4, rowCnt).getContents().trim());
			Nr_node = Integer.parseInt(rs.getCell(5, rowCnt).getContents().trim());
			Nend = Integer.parseInt(rs.getCell(6, rowCnt).getContents().trim());
			NT = Integer.parseInt(rs.getCell(7, rowCnt).getContents().trim());
			Nroute = Integer.parseInt(rs.getCell(8, rowCnt).getContents().trim());
			rowCnt += 4;

			/*
			 * 子系统管段数据表格 Pipe.No 起点号I0 终点号J0 长度LP 直径DP 摩阻系数 起端标高 终端标高 1 0 1 28.5
			 * 0.3 0.017 3.894 3.842 2 1 2 32 0.3 0.017 3.842 3.784 3 2 3 28.6
			 * 0.3 0.017 3.784 3.733 4 3 4 25.4 0.3 0.017 3.733 3.687 5 4 5 24.7
			 * 0.3 0.017 3.687 3.643 6 5 6 23.5 0.3 0.017 3.643 3.601 7 6 7 30.4
			 * 0.3 0.017 3.601 3.546 8 8 7 15.5 0.3 0.017 3.731 3.171 9 9 6 4.3
			 * 0.3 0.017 3.886 3.7
			 */
			I0 = new int[NP];
			J0 = new int[NP];
			lp = new double[NP];
			dpl = new double[NP];
			slp = new double[NP];
			ZJup = new double[NP];
			ZJdw = new double[NP];
			for (int j = 0; j < NP; j++)
			{
				I0[j] = Integer.parseInt(rs.getCell(1, rowCnt + j).getContents().trim());
				J0[j] = Integer.parseInt(rs.getCell(2, rowCnt + j).getContents().trim());
				lp[j] = Double.parseDouble(rs.getCell(3, rowCnt + j).getContents().trim());
				dpl[j] = Double.parseDouble(rs.getCell(4, rowCnt + j).getContents().trim());
				slp[j] = Double.parseDouble(rs.getCell(5, rowCnt + j).getContents().trim());
				ZJup[j] = Double.parseDouble(rs.getCell(6, rowCnt + j).getContents().trim());
				ZJdw[j] = Double.parseDouble(rs.getCell(7, rowCnt + j).getContents().trim());
			}
			rowCnt += NP;
			rowCnt += 3;

			/*
			 * 子系统节点数据表格节点No 汇水面积ha 径流系数 地面标高 井底标高 1 3.5 0.6 5.244 暂未用到 2 3.5
			 * 0.6 5.191 3 3.5 0.6 5.177 4 3.5 0.6 5.208 5 3.5 0.6 5.221 6 3.5
			 * 0.6 5.201 7 3.5 0.6 5.2 8 3.5 0.6 5.121 9 3.5 0.6 5.131 10 3.5
			 * 0.6 5.186
			 */
			Aj = new double[NN];
			Acoef = new double[NN];
			Hj = new double[NN];
			for (int j = 0; j < NN; j++)
			{
				Aj[j] = Double.parseDouble(rs.getCell(1, rowCnt + j).getContents().trim());
				Acoef[j] = Double.parseDouble(rs.getCell(2, rowCnt + j).getContents().trim());
				Hj[j] = Double.parseDouble(rs.getCell(3, rowCnt + j).getContents().trim());
			}
			rowCnt += NN;
			rowCnt += 3;

			/*
			 * 管网路径数&路径节点号节点序号 1 2 3 4 5 6 7 8 1 0 1 2 3 4 5 6 7 2 8 7 -99 -99
			 * -99 -99 -99 -99 3 9 6 -99 -99 -99 -99 -99 -99
			 */
			Mroute = new int[Nstart][Nr_node];
			for (int j = 0; j < Nstart; j++)
			{
				for (int k = 0; k < Nr_node; k++)
				{
					Mroute[j][k] = Integer.parseInt(rs.getCell(k + 1, rowCnt + j).getContents().trim());
				}
			}
			rowCnt += Nstart;
			rowCnt += 3;

			/*
			 * 子系统分支路径管段数据矩阵 倒序pipe branches-reverse order 节点序号 1 2 3 4 5 6 7 1
			 * 6 5 4 3 2 1 0 2 7 -99 -99 -99 -99 -99 -99 3 8 -99 -99 -99 -99 -99
			 * -99
			 */
			Mbranch = new int[Nstart][Npline];
			for (int j = 0; j < Nstart; j++)
			{
				for (int k = 0; k < Npline; k++)
				{
					Mbranch[j][k] = Integer.parseInt(rs.getCell(k + 1, rowCnt + j).getContents().trim());
				}
			}
			// ----临界水深计算变量----
			double sita0 = 3.0, eps = 0.001, alfa = 0.5;
			double Ad0, qkpmax, Hwdwkp, yykp, sita, cons_b, sita_s = 0, sita_c, fsita, dfdsita, dfsita, ssita = 0, csita = 0, hyd_A, hafsita, shafsita = 0, chafsita, sita_p = 0;
			// 中间变量
			int i, j, k, ik, jk, it, k1, kp, in1, in2, in3, NR1, NR2, ii, Nprt, iprt1, iprt2;
			double H00, ycd0;
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
			double[][] Hwj = new double[NT][NN];

			// ----------------------------------------------------------------------------------------------------------
			String FileName = "";
			if (gjId != null)
			{
				FileName = gjId.substring(0, 12) + ".txt";
			}
			else
			{
				FileName = subSys + ".txt";
			}
			String FilePath = "./www/DPP-LOCAL/DPP-LOCAL-WEB/files/analogValue/";
			FileOutputStream fs = new FileOutputStream(new File(FilePath + FileName));
			PrintStream printStream = new PrintStream(fs);
			printStream.println(FileName);

			DecimalFormat df = new DecimalFormat("##.####");
			DecimalFormat df1 = new DecimalFormat("######.##");
			// --输出数据文件开始---
			// ================= 赋初值 ===============================
			for (i = 0; i < NT; i++)
			{
				for (j = 0; j < NN; j++)
					sumAj[i][j] = 0;
			}
			for (i = 0; i < NT; i++)
			{
				for (j = 0; j < NN; j++)
					sumqj[i][j] = 0;
			}
			for (i = 0; i < NN; i++)
			{
				for (j = 0; j < NN; j++)
				{
					if (i == j)
					{
						Tnode[i][j] = 0;
					}
					else
					{
						Tnode[i][j] = -99;
					}
				}
			}
			for (i = 0; i < NN; i++)
			{
				for (j = 0; j < NN; j++)
					sumTnode[i][j] = 0;
			}
			// ==================Tnode-sumTnode=========================
			for (i = 0; i < NP; i++)
				vp[i] = vp0;
			for (kp = 0; kp < NP; kp++)
			{
				in1 = I0[kp];
				in2 = J0[kp];
				Tnode[in1][in2] = lp[kp] / vp[kp] / 60;
				slop[kp] = (ZJup[kp] - ZJdw[kp]) / lp[kp];
			}
			for (i = 0; i < Nroute; i++)
			{
				for (j = 0; j < Nr_node; j++)
				{
					in1 = Mroute[i][j];
					if (in1 >= 0)
					{
						for (k = j + 1; k < Nr_node; k++)
						{
							in2 = Mroute[i][k - 1];
							in3 = Mroute[i][k];
							if (in3 >= 0)
							{
								sumTnode[in1][in3] = sumTnode[in1][in2] + Tnode[in2][in3];
							}
						}
					}
				}
			}
			// System.out.println("pipe no.  I0    J0");
			for (i = 0; i < NP; i++)
			{
				// System.out.printf("%6d%6d%6d", i, I0[i], J0[i]);
				// System.out.println();
			}
			printStream.println();
			printStream.print(" ip=");
			for (i = 0; i < NP; i++)
			{
				printStream.printf("%4d", i);
			}
			printStream.println();
			printStream.print(" I0=");
			for (i = 0; i < NP; i++)
			{
				printStream.printf("%4d", I0[i]);
			}
			printStream.println();
			printStream.print(" J0=");
			for (i = 0; i < NP; i++)
			{
				printStream.printf("%4d", J0[i]);
			}
			printStream.println();
			printStream.println();
			printStream.println("===========  print Mroute[i][j]");
			for (i = 0; i < Nroute; i++)
			{
				for (j = 0; j < Nr_node; j++)
				{
					printStream.printf("%6d", Mroute[i][j]);
				}
				printStream.println();
			}
			printStream.println();
			printStream.println("===========  print Mbranch[i][j]");
			for (i = 0; i < Nstart; i++)
			{
				for (j = 0; j < Npline; j++)
				{
					printStream.printf("%6d", Mbranch[i][j]);
				}
				printStream.println();
			}
			printStream.println("===========  print Tnode[i][j]");
			printStream.println("====j=  ");
			printStream.println("      ");
			for (j = 0; j < NN; j++)
			{
				printStream.printf("%6d", j);
			}
			printStream.println();
			for (i = 0; i < NN; i++)
			{
				if (i < 10)
				{
					printStream.print("i=" + i + "   ");
				}
				else
				{
					printStream.print("i=" + i + "  ");
				}
				for (j = 0; j < NN; j++)
				{
					if (Tnode[i][j] < 0.0)
					{
						printStream.print("      ");
					}
					else
					{
						printStream.printf("%6.2f", Tnode[i][j]);
					}
				}
				printStream.println();
			}
			printStream.println();
			printStream.println("===========  print sumTnode[i][j]");
			printStream.print("==j=  ");
			for (j = 0; j < NN; j++)
			{
				printStream.printf("%6d", j);
			}
			printStream.println();
			for (i = 0; i < NN; i++)
			{
				printStream.print("i=" + i + "   ");
				for (j = 0; j < NN; j++)
				{
					if (sumTnode[i][j] <= 0.0)
					{
						printStream.print("      ");
					}
					else
					{
						printStream.printf("%6.2f", sumTnode[i][j]);
					}
				}
				printStream.println();
			}
			// ================= 管网准稳态流动模拟============================
			// -------------------动态模拟流量计算-----------------------------
			// ----------------节点汇水面积(ha)和汇水流量(m3/sec)计算--------
			printStream.println();
			printStream.println("===========  管网动态模拟计算      重现期＝ " + P_simu + "  年   时段数＝ " + NT + "       终点水位＝ " + Hw_end + "  m  =========");
			// 芝加哥过程线--rainfall intensity at every time step--
			AA = A1 + A1 * C_storm * Math.log(P_simu) / 2.303;
			for (it = 0; it < NT; it++)
			{
				if (it <= NR)
				{
					dtnt = dt * (float) (it);
					tbb = dt * (float) (NR) - dtnt;
					XX1 = AA * ((1.0 - n_storm) * tbb / rc + b_storm);
					XX2 = Math.pow((tbb / rc + b_storm), (n_storm + 1.0));
				}
				else
				{
					dtnt = dt * (float) (it);
					taa = dtnt - dt * (float) (NR);
					XX1 = AA * ((1.0 - n_storm) * taa / (1.0 - rc) + b_storm);
					XX2 = Math.pow((taa / (1.0 - rc) + b_storm), (n_storm + 1.0));
				}
				XX[it] = XX1 / XX2;
				qit[it] = 167.0 * XX[it] / 1000.0;
			}
			NR1 = NR - 1;
			NR2 = NR + 1;
			qit[NR] = (qit[NR] + qit[NR - 1] + qit[NR + 1]) / 3.0;
			printStream.println();
			printStream.println("    it      dtnt      XX[it]     qit[it]");
			for (it = 0; it < NT; it++)
			{
				dtnt = dt * (float) (it);
				printStream.printf("%6d%10.2f%12.6f%12.6f", it, dtnt, XX[it], qit[it]);
				printStream.println();
			}
			printStream.println();
			for (it = 0; it < NT; it++)
			{
				dtnt = dt + dt * (float) (it);
				for (j = 0; j < NN; j++)
				{
					sumAj[it][j] = Aj[j];
					sumqj[it][j] = Aj[j] * qit[it] * Acoef[j];
					for (i = 0; i < NN; i++)
					{
						if (sumTnode[i][j] > 0 && sumTnode[i][j] < dtnt)
						{
							sumAj[it][j] = sumAj[it][j] + Aj[i];
							sumqj[it][j] = sumqj[it][j] + Aj[i] * qit[it] * Acoef[i];
						}
					}
				}
			}
			printStream.println("  sumAj[it][j]=");
			for (it = 0; it < NT; it++)
			{
				for (j = 0; j < NN; j++)
				{
					printStream.printf("%8.2f", sumAj[it][j]);
				}
				printStream.println();
			}
			printStream.println();
			printStream.println("  sumqj[it][j]=");
			for (it = 0; it < NT; it++)
			{
				for (j = 0; j < NN; j++)
				{
					printStream.printf("%8.2f", sumqj[it][j]);
				}
				printStream.println();
			}
			printStream.println();
			// ---------------------------------------------------------------
			for (it = 0; it < NT; it++)
			{
				for (i = 0; i < NN; i++)
				{
					overflow[it][i] = 0.0;
					Hw_over[it][i] = 0.0;
				}
			}
			for (it = 0; it < NT; it++)
			{
				for (j = 0; j < NP; j++)
				{
					qpt[it][j] = -99.0;
					qqkp[it][j] = 0.0;
				}
			}
			// ---------------------------------------------------------------
			for (it = 0; it < NT; it++)
			{
				printStream.print(" it=" + it + "  qpt[it][k]=");
				for (j = 0; j < NN; j++)
				{
					for (k = 0; k < NP; k++)
					{
						if (I0[k] == j)
						{
							qpt[it][k] = sumqj[it][j];
							printStream.printf("%8.2f", qpt[it][k]);
						}
					}
				}
				printStream.println();
				for (ik = 0; ik < Nstart; ik++)
				{
					for (jk = 0; jk < Npline; jk++)
					{
						kp = Mbranch[ik][jk];
						if (kp >= 0)
						{
							if (J0[kp] == Nend)
							{
								Hwdw[it][kp] = Hw_end;
								if (1 == Iprt)
								{
									printStream.println("   it= " + it + "   kp= " + kp + "  Hwdm= " + Hwdw[it][kp] + "  Hw_end= " + Hw_end);
								}
							}
							else
							{
								for (k1 = 0; k1 < NP; k1++)
								{
									if (I0[k1] == J0[kp]) Hwdw[it][kp] = Hwup[it][k1];
								}
							}
							Ad0 = 0.7854 * Math.pow(dpl[kp], 2.0);
							hdj0 = ZJdw[kp] + dpl[kp];
							if (Hwdw[it][kp] >= hdj0)
							{
								if (1 == Iprt)
								{
									printStream.println("   it= " + it + "   kp= " + kp + "  Hwdm= " + df.format(Hwdw[it][kp]) + "  淹没出流 ");
								}
								hdcc0[it][kp] = 1.0;
								rid[it][kp] = dpl[kp] / 4.0;
								vpt[it][kp] = qpt[it][kp] / Ad0;
								slopt[it][kp] = 10.29 * Math.pow(slp[kp], 2.0) * Math.pow(qpt[it][kp], 2.0) / Math.pow(dpl[kp], 5.333);
								Hwup[it][kp] = Hwdw[it][kp] + slopt[it][kp] * lp[kp];
								if (Hwup[it][kp] >= Hj[I0[kp]])
								{
									Hwup[it][kp] = Hj[I0[kp]];
									slopt[it][kp] = (Hwup[it][kp] - Hwdw[it][kp]) / lp[kp];
									if (slopt[it][kp] < 0.0)
									{
										slopt[it][kp] = Math.abs(slopt[it][kp]);
									}
									vpt[it][kp] = Math.pow(rid[it][kp], 0.6667) * Math.pow(slopt[it][kp], 0.5) / slp[kp];
									qqkp[it][kp] = vpt[it][kp] * Ad0;
									if (qqkp[it][kp] < 0.0)
									{
										qqkp[it][kp] = Math.abs(qqkp[it][kp]);
									}
								}
							}
							else
							{
								if (Iprt == 1)
								{
									printStream.println("   it= " + it + "   kp= " + kp + "  Hwdw= " + Hwdw[it][kp] + "  非淹没出流 ");
								}
								// --20161018修改开始---采用临界水深简化算法-----------------------
								qkpmax = 2.699 * Math.pow(dpl[kp], 2.5);
								if (qpt[it][kp] > qkpmax)
								{
									if (Iprt == 1)
									{
										printStream.println("   it= " + it + "   kp= " + kp + "  qkpmax= " + qkpmax + "  非淹没满管出流 ");
									}
									vpt[it][kp] = qpt[it][kp] / Ad0;
									H00 = Math.pow(vpt[it][kp], 2.0) / 13.72;
									Hwdw[it][kp] = ZJdw[kp] + dpl[kp] + H00;
									hdcc0[it][kp] = 1.0;
									rid[it][kp] = dpl[kp] / 4.0;
									slopt[it][kp] = 10.29 * Math.pow(slp[kp], 2.0) * Math.pow(qpt[it][kp], 2.0) / Math.pow(dpl[kp], 5.333);
									Hwup[it][kp] = Hwdw[it][kp] + slopt[it][kp] * lp[kp];
									if (Hwup[it][kp] >= Hj[I0[kp]])
									{
										Hwup[it][kp] = Hj[I0[kp]];
										slopt[it][kp] = (Hwup[it][kp] - Hwdw[it][kp]) / lp[kp];
										if (slopt[it][kp] < 0.0)
										{
											slopt[it][kp] = Math.abs(slopt[it][kp]);
										}
										vpt[it][kp] = Math.pow(rid[it][kp], 0.6667) * Math.pow(slopt[it][kp], 0.5) / slp[kp];
										qqkp[it][kp] = vpt[it][kp] * Ad0;
										if (qqkp[it][kp] < 0.0)
										{
											qqkp[it][kp] = Math.abs(qqkp[it][kp]);
										}
									}
								}
								else
								{
									if (Iprt == 1)
									{
										printStream.println("   it= " + it + "   kp= " + kp + "  Hwdm= " + Hwdw[it][kp] + "  非淹没非满管出流 ");
									}
									// ==20161018修改开始---采用临界水深简化公式--------zhou-p21------
									ycd0 = qpt[it][kp] / 2.983 / Math.pow(dpl[kp], 2.5);
									hdcc0[it][kp] = Math.pow(ycd0, 0.513);
									sita = 2.0 * Math.acos(1.0 - 2.0 * hdcc0[it][kp]);
									rid[it][kp] = 0.25 * dpl[kp] * (sita - Math.sin(sita)) / sita;
									vpt[it][kp] = Math.pow(rid[it][kp], 0.6667) * Math.pow(slop[kp], 0.5) / slp[kp];
								}
								// ---for(k=0;k<N;k++)结束---20160907修改结束---临界水深算法--------------
								Hwdwkp = ZJdw[kp] + hdcc0[it][kp] * dpl[kp];
								if (Hwdwkp >= Hwdw[it][kp])
								{
									Hwdw[it][kp] = Hwdwkp;
								}
								else
								{
									yykp = Hwdw[it][kp] - ZJdw[kp];
									if (yykp > dpl[kp])
									{
										yykp = dpl[kp];
									}
									sita = 2.0 * Math.acos(1.0 - 2.0 * yykp / dpl[kp]);
									hdcc0[it][kp] = yykp / dpl[kp];
									rid[it][kp] = 0.25 * dpl[kp] * (sita - Math.sin(sita)) / sita;
									vpt[it][kp] = Math.pow(rid[it][kp], 0.6667) * Math.pow(slop[kp], 0.5) / slp[kp];
								}
								Hwup[it][kp] = Hwdw[it][kp] + slop[kp] * lp[kp];
							}
							// ------- 输出it计算结果 ----------
							if (Iprt == 1)
							{
								printStream.println("   it= " + it + "   kp= " + kp + "   I0[kp]= " + I0[kp] + "  Hwdm= " + Hwdw[it][kp] + "  Hwup= " + Hwup[it][kp] + "  Hj= " + Hj[I0[kp]] + "  hdcc0= " + hdcc0[it][kp] + "  qpt= " + qpt[it][kp] + "  qqkp= " + qqkp[it][kp] + "  vpt= " + vpt[it][kp]);
							}
						}
					}
				}
				printStream.println();

				printStream.println("    it   管段号  I0   J0 管径dpl     管段qp   水力半径R  充满度 流速(m/s)  上游水位  下游水位  上管底高  下管底高  管段坡度  上地面高");
				for (i = 0; i < NP; i++)
				{
					printStream.printf("%6d%6d%6d%5d%8.2f%12.3f%10.3f%8.3f%10.3f%10.3f%10.3f%10.3f%10.3f%10.5f%10.3f", it, i, I0[i], J0[i], dpl[i], qpt[it][i], rid[it][i], hdcc0[it][i], vpt[it][i], Hwup[it][i], Hwdw[it][i], ZJup[i], ZJdw[i], slop[i], Hj[I0[i]]);
					printStream.println();
				}
				printStream.println();
				// -------------- 开始计算溢流节点 ---------------
				for (i = 0; i < NP; i++)
				{
					k = J0[i];
					if (k == Nend)
					{
						Hwj[it][k] = Hwdw[it][i];
					}
					{
						j = I0[i];
						Hwj[it][j] = Hwup[it][i];
						if (Hwup[it][i] == Hj[j])
						{
							overflow[it][j] = overflow[it][j] + (qpt[it][i] - qqkp[it][i]) * dt * 60.0;
							Hw_over[it][j] = csf * overflow[it][j] / Aj[j] / 10000.0 * 1000.0;

						}
						if (Hwup[it][i] < Hj[j] && overflow[it][j] > 0.0)
						{
							overflow[it][j] = overflow[it - 1][j] * 0.90;
							Hw_over[it][j] = csf * overflow[it][j] / Aj[j] / 10000.0 * 1000.0;
						}
					}
					if (it > NR && Hw_over[it][j] <= 5.0)
					{
						overflow[it][j] = 0.0;
						Hw_over[it][j] = 0.0;
					}
				}
				// ------------------ 计算溢流节点结束 ---------------
			}
			// ----------------屏幕输出计算结束------
			//System.out.println("------ 模型计算全部完成 ------");
			// ---------------------- 输出管段充满度计算结果 ---------------
			printStream.println(" ======== 时段管段充满度 ========");
			Nprt = NP / Nprtc + 1;
			for (ii = 0; ii < Nprt; ii++)
			{
				{
					iprt1 = ii * Nprtc;
					iprt2 = iprt1 + Nprtc;
					if (iprt2 > NP)
					{
						iprt2 = NP;
					}
				}
				printStream.print("  i=    ");
				for (i = iprt1; i < iprt2; i++)
				{
					if (i < 10)
					{
						printStream.print("    " + i + "   ");
					}
					else
					{
						printStream.print("   " + i + "   ");
					}
				}
				printStream.println();
				printStream.println();
				for (it = 0; it < NT; it++)
				{
					if (it < 10)
					{
						printStream.print(" " + it + "   ");
					}
					else
					{
						printStream.print(it + "   ");
					}
					for (i = iprt1; i < iprt2; i++)
					{
						printStream.printf("%8.3f", hdcc0[it][i]);
					}
					printStream.println();
				}
			}
			//
			// ------------------- 输出节点水位计算结果 ---------------
			printStream.println(" ======== 时段节点水位 ========");
			Nprt = NN / Nprtc + 1;
			for (ii = 0; ii < Nprt; ii++)
			{
				{
					iprt1 = ii * Nprtc;
					iprt2 = iprt1 + Nprtc;
					if (iprt2 > NN)
					{
						iprt2 = NN;
					}
				}
				printStream.print("  i=    ");
				for (i = iprt1; i < iprt2; i++)
				{
					if (i < 10)
					{
						printStream.print("    " + i + "   ");
					}
					else
					{
						printStream.print("   " + i + "  ");
					}
				}
				printStream.println();
				printStream.println("it=");
				for (it = 0; it < NT; it++)
				{
					if (it < 10)
					{
						printStream.print(" " + it + "   ");
					}
					else
					{
						printStream.print(it + "   ");
					}
				}
			}
			
			//***********组织数据，传到页面用于显示********
			for (it = 0; it < NT; it++)
			{
				String WaterLevNew = "";
				for (i = 0; i < NN; i++)
				{
					if (gjId != null && i == SubgjId)
					{
						WaterAccGj += df1.format(Hwj[it][i]) + "|";
					}
					WaterLevNew += df1.format(Hwj[it][i]) + "|";
				}
				WaterLev[it] = WaterLevNew;
			}
			//*************************************
			// ------------------ 输出节点溢流计算结果 ---------------
			printStream.println(" ======== 时段节点积水量(m3) ========");
			Nprt = NN / Nprtc + 1;
			for (ii = 0; ii < Nprt; ii++)
			{
				{
					iprt1 = ii * Nprtc;
					iprt2 = iprt1 + Nprtc;
					if (iprt2 > NN)
					{
						iprt2 = NN;
					}
				}
				printStream.print("  i=    ");
				for (i = iprt1; i < iprt2; i++)
				{
					if (i < 10)
					{
						printStream.print("    " + i + "   ");
					}
					else
					{
						printStream.print("   " + i + "  ");
					}
				}
				printStream.println();
				printStream.println("it=");
				for (it = 0; it < NT; it++)
				{
					if (it < 10)
					{
						printStream.println(" " + it + "   ");
					}
					else
					{
						printStream.println(it + "   ");
					}
					for (i = iprt1; i < iprt2; i++)
					{
						if (overflow[it][i] <= 0.0)
						{
							printStream.print("        ");
						}
						else
						{
							printStream.printf("%8.2f", overflow[it][i]);
						}
					}
					printStream.println();
					String WaterAccNew = "";
					for (i = 0; i < NN; i++)
					{
						if (overflow[it][i] <= 0.0)
						{
							WaterAccNew += 0 + "|";
						}
						else
						{
							WaterAccNew += df1.format(overflow[it][i]) + "|";
						}
					}
					WaterAcc[it] = WaterAccNew;
				}
			}
			
			//***********组织数据，传到页面用于显示********
			for (it = 0; it < NT; it++)
			{
				String WaterAccNew = "";
				for (i = 0; i < NN; i++)
				{
					if (overflow[it][i] <= 0.0)
					{
						WaterAccNew += 0 + "|";
					}
					else
					{
						WaterAccNew += df1.format(overflow[it][i]) + "|";
					}
				}
				WaterAcc[it] = WaterAccNew;
			}
			//*********************************************
			printStream.println(" ======== 时段节点积水深度(mm) ========");
			Nprt = NN / Nprtc + 1;
			for (ii = 0; ii < Nprt; ii++)
			{
				{
					iprt1 = ii * Nprtc;
					iprt2 = iprt1 + Nprtc;
					if (iprt2 > NN)
					{
						iprt2 = NN;
					}
				}
				printStream.print("  i=    ");
				for (i = iprt1; i < iprt2; i++)
				{
					if (i < 10)
					{
						printStream.print("    " + i + "   ");
					}
					else
					{
						printStream.print("   " + i + "   ");
					}
				}
				printStream.println();
				printStream.println("it=");
				for (it = 0; it < NT; it++)
				{
					if (it < 10)
					{
						printStream.print(" " + i + "   ");
					}
					else
					{
						printStream.print(i + "   ");
					}
					for (i = iprt1; i < iprt2; i++)
					{
						if (overflow[it][i] <= 0.0)
						{
							printStream.print("        ");
						}
						else
						{
							printStream.printf("%8.2f", Hw_over[it][i]);
						}
					}
					printStream.println();
				}
			}
		}
		catch (Exception e)
		{
			e.printStackTrace();
		}
		if (AnalogWaterType.equals("WaterAccGj"))
		{
			return WaterAccGj;
		}
		else if (AnalogWaterType.equals("WaterAcc"))
		{
			String WaterAccList = "";
			for (int i = 0; i < WaterAcc.length; i++)
			{
				WaterAccList += subSys.substring(7, 12) + WaterAcc[i] + ";";
			}
			return WaterAccList;
		}
		else if (AnalogWaterType.equals("WaterLev"))
		{
			String WaterLevList = "";
			for (int i = 0; i < WaterLev.length; i++)
			{
				WaterLevList += subSys.substring(7, 12) + WaterLev[i] + ";";
			}
			return WaterLevList;
		}
		return "";
	}

	private String		FileSaveRoute;
	private String		File_Name;
	private String		Sid;

	private String		AnalogWaterType;
	private String[]	WaterAcc	= new String[60];
	private String[]	WaterLev	= new String[60];
	private String		WaterAccGj	= "";
}
