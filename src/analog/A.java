package analog;

public class A {

//20160825---- 多子系统雨水排水管网
//20150808new雨水管网过程模拟计算程序-芝加哥过程线
//－－SQ LIU, TONGJI UNIVERSITY, 25 AUG 2016 
//

//
public static void main()
{

//管网基础数据：子系统数，最大管段数，最大节点数，最大管道起点数，路径最大管段数，最大计算次数，中间结果写文件指针
int Nbr=8,NP=116,NN=117,Nstart=10,Npline=10,Nmax=20,Iprt=0; 
//模拟时段数，芝加哥峰点时段位置（r=0.375），管道路径数，路径最大节点数，终点节点号
int NT=60,NR=23,Nroute=10,Nr_node=10,Nend=8;
int ip,it,k1,kp,in1,in2,in3;
double dtnt,taa,tbb,AA,XX1,XX2;
double []XX = new double[NT];
double []qit = new double[NT];
double [][]sumqj = new double[NT][NN];
double [][]sumAj = new double[NT][NN];
double [][]Tnode = new double[NN][NN];
double [][]sumTnode = new double[NN][NN];
int i,j,k,ik,jk,k00,kkk1;
double hdc_min,hdc_max;
double []hdcc = new double[2];
double []fun_hd = new double[2];
double fun_hd0;
double [][]hdcc0= new double[NT][NP];
double hda,hdb,hdc = 0,hde = 0,hdf,hdd;
//管段流速（m/s）
double []vp = new double[NP];
double []slop = new double[NP];
double Ad0,AD,RD,vp0=1.0;
double [][]qpt = new double[NT][NP];
double [][]vpt=new double[NT][NP];
double [][]rid = new double[NT][NP];
double [][]slopt=new double[NT][NP];
double[][] Hwup = new double[NT][NP];
double [][]Hwdw =new double[NT][NP];
double [][]qjt = new double[NT][NN];
double [][]overflow = new double[NT][NN];
double[][] Hw_over=new double[NT][NN];
double hdj0,ARD23;
double []totalflow=new double[NN];
double[] totalHw =new double[NN];
//暴雨公式shanghai storm water formular: (A1+C*lgP)/(t+b)**n---ln(N)=2.303log(N)---出口水位（m）
double A1=17.53, C_storm=0.95, tmin=10, b_storm=11.77, P_simu=10,n_storm=0.88, dt=2.0,rc=0.375,Hw_end=0;
//节点汇水面积(ha)
double []Aj = new double[NN];
//节点汇水面积径流系数
double []Acoef = new double[NN];
//节点地面标高（m） [NN=23]
double[] Hj = new double[NN];
//管网起始节点号和起始节点管底埋深<m>
//int NJstart [Nstart]={ };
//double HJstart [Nstart]={ };
//管网路径数和路径节点号(-99表示空节点)
int [][]Mroute=new int [Nroute][Nr_node];
//pipe branches-reverse order
int [][]Mbranch=new int[Nstart][Npline];
//管段上游节点号I0,下游节点号J0，管段长度(m),摩|阻系数
int []I0 =new int[NP];
int []J0 = new int[NP];
double []lp =new double[NP];
//double slp [NP]={};
double []slp = new double[NP];
//管段直径(m)，上游管底高程(m)，下游管底高程(m)
double []dpl = new double[NP];
double []ZJup = new double[NP];
double []ZJdw = new double[NP];
//----------------------------------------------------------------------------------------------------------
//
//--输出数据文件开始---
//ofstream outfile;
//outfile.open("20160808-雨水管网过程模拟-芝加哥过程线-10-60.txt");
System.out.println("20160808-雨水管网过程模拟-芝加哥过程线-10-60.txt");
//================= 赋初值 ===============
//for(=0;i<Nbr;++){
	for(i=0;i<NT;i++){for(j=0;j<NN;j++) sumAj [i][j]=0;}
	  for(i=0;i<NT;i++){for(j=0;j<NN;j++) sumqj [i][j]=0;}
	  for(i=0;i<NN;i++){for(j=0;j<NN;j++) {if(i==j){Tnode [i][j]=0;}else{Tnode [i][j]=-99;}}}
	  for(i=0;i<NN;i++){for(j=0;j<NN;j++) sumTnode [i][j]=0;}
		//}
//==================Tnode-sumTnode=========================
//for(=0;i<Nbr;++){
for (i=0;i<NP;i++) vp [i]=vp0; 
 for(kp=0;kp<NP;kp++)
	   {in1=I0 [kp];
      in2=J0 [kp];
		Tnode [in1][in2]=lp [kp]/vp [kp]/60; 
		slop [kp]=(ZJup [kp]-ZJdw [kp])/lp [kp];
		}
//
	 for(i=0;i<Nroute;i++) 
	 {for(j=0;j<Nr_node;j++) 
	 {in1=Mroute [i][j];
    if(in1>=0)
	 {for(k=j+1;k<Nr_node;k++) 
	  { in2=Mroute [i][k-1]; in3=Mroute [i][k]; 
     if(in3>=0)
	   {sumTnode [in1][in3]=sumTnode [in1][in2]+Tnode [in2][in3];
	 }}}}}
//}  
//
//=====print Mroute[i][j], Tnode, sumTnode,Mbranch[i][j]====
 //cout<<"pipe no.  I0    J0"<<endl;
	 System.out.println("pipe no.  I0    J0");
 //for(=0;i<Nbr;++){
for(i=0;i<NP;i++)
  {//s.Format("%6d%6d%6d\n", i,I0 [i],J0 [i]);
	  // cout<<s;
	System.out.print(i+I0 [i]+J0 [i]);
	}
   //outfile<<endl; 
System.out.println();
//
 	// outfile<<" ip="; 
System.out.print(" ip=");
	 for(i=0;i<NP;i++)
	  { //s.Format("%4d",i);
	   //outfile<<s;
		 System.out.print(i);
	  }
  // outfile<<endl; 
	 System.out.println();
//
 	// outfile<<" I0="; 
	 System.out.print(" I0=");
	 for(i=0;i<NP;i++)
	  { //s.Format("%4d",I0 [i]);
	   //outfile<<s;
		 System.out.print(I0 [i]);
	  }
  // outfile<<endl; 
	 System.out.println();
//
 	 //outfile<<" J0="; 
 	 System.out.print(" J0=");
	 for(i=0;i<NP;i++)
	  { //s.Format("%4d",J0 [i]);
	   //outfile<<s;
		 System.out.print(J0 [i]);
	  }
   //outfile<<endl; 
	 System.out.println();
//}//
//
   //outfile<<endl; 
   //outfile<<"===========  print Mroute [i][j]"<<endl; 
	 System.out.println("===========  print Mroute [i][j]");
	//for(=0;i<Nbr;++) {
for(i=0;i<Nroute;i++) 
	 {for(j=0;j<Nr_node;j++) 
	 {	 //s.Format("%6d",Mroute [i][j]);
	  // outfile<<s;
		 System.out.print(Mroute [i][j]);
	 }
	 //outfile<<endl; 
	 System.out.println();
	 }
//}//
//
   //outfile<<endl; 
  // outfile<<"===========  print Mbranch [i][j]"<<endl; 
	System.out.println("===========  print Mbranch [i][j]");
	// for(=0;i<Nbr;++) {
for(i=0;i<Nstart;i++) 
	 {for(j=0;j<Npline;j++) 
	 {	 //s.Format("%6d",Mbranch [i][j]);
	   //outfile<<s;
		 System.out.print(Mbranch [i][j]);
	 }
	// outfile<<endl; 
	 System.out.println();
	 }
//}//
//
  // outfile<<"===========  print Tnode [i][j]"<<endl; 
  // outfile<<"====j=  "<<endl; 
  // outfile<<"      "; 
	System.out.println("===========  print Tnode [i][j]");
	System.out.println("====j=  ");
	//for(=0;i<Nbr;++) {
for(j=0;j<NN;j++) 
	 {	 //s.Format("%6d",j);
	  // outfile<<s;
	System.out.print(j);
	 }
      //outfile<<endl; 
System.out.println();
//
   for(i=0;i<NN;i++)
	 { 
	   if(i<10){
		   //outfile<<"i="<<i<<"   ";
		   System.out.print("i=" + i + "   ");
		   }else{
			   //outfile<<"i="<<i<<"  ";
			   System.out.print("i=" + i + "   ");
	 }
	 for(j=0;j<NN;j++){
		 if(Tnode [i][j]<0.0) {
			 //outfile<<"      ";
			 } else{	
				 //s.Format("%6.2f",Tnode [i][j]);
				 //outfile<<s;
				 System.out.print(Tnode [i][j]);
	 }}
      //outfile<<endl;
	 System.out.println();
	 }
//}//
//
   //outfile<<endl; 
   //outfile<<"===========  print sumTnode [i][j]"<<endl; 
  // outfile<<"==j=  "; 
//   outfile<<"      "; 
   System.out.println("===========  print sumTnode [i][j]");
   System.out.println("==j=       ");
//	for(=0;i<Nbr;++){
for(j=0;j<NN;j++) 
	 {	 //s.Format("%6d",j);
	  // outfile<<s;
	System.out.print(j);
	 }
      //outfile<<endl;
System.out.println();
//
   for(i=0;i<NN;i++)
	 { //outfile<<"i="<<i<<"   ";
	   
	 System.out.println("i=" + i + "   ");
	 for(j=0;j<NN;j++)
	 {if(sumTnode [i][j]<=0.0) {
		 //outfile<<"      ";
		 System.out.print("  ");
		 } else
	 
	 {	 //s.Format("%6.2f",sumTnode [i][j]);
	   //outfile<<s;
			 System.out.print(sumTnode [i][j]);
	 }}
     // outfile<<endl; 
	 System.out.println();
	 }
//}//
//================= 管网稳态流动模拟===================
//
//-----------动态模拟流量计算------------------------------------
//----------------节点汇水面积(ha)和汇水流量(m3/sec)计算--------
	//outfile<<endl; 
	//outfile<<"===========  管网动态模拟计算      重现期＝ "<<P_simu<<"  年   时段数＝ "<<NT<<"       终点水位＝ "<<Hw_end<<"  m  ========="<<endl;
System.out.println("\n===========  管网动态模拟计算      重现期＝ "+P_simu+"  年   时段数＝ "+NT+"       终点水位＝ "+Hw_end+"  m  =========");
//
//----------rainfall intensity at every time step-----------
//xxxxxxx
//芝加哥过程线
 AA=A1+A1*C_storm*Math.log(P_simu)/2.303;
 for (it=0;it<NT;it++)
 {if (it<=NR) 
 {dtnt=dt * (float)(it);
 tbb=dt* (float)(NR) - dtnt;
 XX1=AA*((1.0-n_storm)*tbb/rc+b_storm);
 XX2=Math.pow((tbb/rc+b_storm),(n_storm+1.0));
 } else
 {dtnt=dt * (float)(it);
 taa= dtnt - dt*(float)(NR);
 XX1=AA*((1.0-n_storm)*taa/(1.0-rc)+b_storm);
 XX2=Math.pow((taa/(1.0-rc)+b_storm),(n_storm+1.0));
 }
 XX[it]=XX1/XX2;
 qit[it]=167.0*XX[it]/1000.0;
 }
//
 //outfile<<endl;
 //outfile<<"    it      dtnt      XX[it]     qit[it]"<<endl;
 System.out.println("\nit\tdtnt\tXX[it]\tqit[it]");
	for(it=0;it<NT;it++)
  {dtnt=dt * (float)(it);
    //s.Format("%6d%10.2lf%12.6lf%12.6lf\n",it,dtnt,XX[it],qit[it]);
	   //outfile<<s;
  System.out.println(it+"\t"+dtnt+"\t"+XX[it]+"\t"+qit[it]);
  }
 // outfile<<endl;
	System.out.println();
//xxxxxxx	
//
//for(=0;i<Nbr;++){
for (it=0;it<NT;it++)
 {dtnt=dt + dt * (float)(it);
 for (j=0;j<NN;j++)
 {sumAj [it][j]=Aj [j]; 
  sumqj [it][j]=Aj [j]*qit[it]*Acoef [j];
 for (i=0;i<NN;i++)
 {if(sumTnode [i][j]>0 && sumTnode [i][j]<dtnt) 
 {sumAj [it][j]=sumAj [it][j]+Aj [i];
  sumqj [it][j]=sumqj [it][j]+Aj [i]*qit[it]*Acoef [i];
 }}}}
//print sumAj[it][j] and sumqj[it][j]
// outfile<<"  sumAj [it][j]="<<endl; 
System.out.println("  sumAj [it][j]=");
 for (it=0;it<NT;it++)
 {for (j=0;j<NN;j++)
 {//s.Format("%8.2lf",sumAj [it][j]);
  //outfile<<s;
	 System.out.print(sumAj [it][j]);
 }
 // outfile<<endl;
 System.out.println();
 }
  //outfile<<endl;
System.out.println("\n  sumqj [it][j]=");
// outfile<<"  sumqj [it][j]="<<endl; 
 for (it=0;it<NT;it++)
 {for (j=0;j<NN;j++)
 {//s.Format("%8.2lf",sumqj [it][j]);
	  // outfile<<s;
	 System.out.print(sumqj [it][j]);
 }
  //outfile<<endl;
 System.out.println();
 }
  //outfile<<endl;
 System.out.println();
//}//
//---------------------------------------------------------------    
 for(it=0;it<NT;it++)  {for(i=0;i<NN;i++) {overflow [it][i]=0.0; Hw_over [it][i]=0.0;}}
 for(it=0;it<NT;it++) { for(j=0;j<NP;j++) qpt [it][j]=-99.0;}
//---------------------------------------------------------------	
// -----000开始
//for(=0;i<Nbr;++){
  for(it=0;it<NT;it++)
//--1
 {//outfile<<" it="<<it<<"  qpt [it][k]=";
	  System.out.print(" it="+it+"  qpt [it][k]=");
  for(j=0;j<NN;j++)
  {for(k=0;k<NP;k++) 
	 {if(I0 [k]==j) 
	  {qpt [it][k]=sumqj [it][j];
    // s.Format("%8.2lf",qpt [it][k]);  outfile<<s;
	  System.out.print(qpt [it][k]);
	  }
   }
	}
  //outfile<<endl;
  System.out.println();
//-------------------???????????????090127?????????????????------------------------
	  for(ik=0; ik<Nstart; ik++) 
//--2
	  {for(jk=0; jk<Npline; jk++)
//--3
	   {kp=Mbranch [ik][jk];
	    if(kp>=0)
//--4
	    {if(J0 [kp]==Nend) {Hwdw [it][kp]=Hw_end;
		//outfile<<"   it= "<<it<<"   kp= "<<kp<<"  Hwdm= "<<Hwdw[it][kp]<<"  Hw_end= "<<Hw_end<<endl; 
	    System.out.println("   it= "+it+"   kp= "+kp+"  Hwdm= "+Hwdw[it][kp]+"  Hw_end= "+Hw_end);
		}
		 else
		 {for(k1=0;k1<NP;k1++) 
		  {if(I0 [k1]==J0 [kp]) Hwdw [it][kp]=Hwup [it][k1];
		 }}
//       
		 Ad0=0.7854*Math.pow(dpl [kp],2.0);
		 hdj0=ZJdw[kp]+dpl [kp];
       if(Hwdw [it][kp]>=hdj0)
		 {hdcc0[it][kp]=1.0;  rid [it][kp]=dpl [kp]/4.0;
	      vpt [it][kp]=qpt[it][kp]/Ad0;
        slopt[it][kp]=10.29*Math.pow(slp [kp],2.0)*Math.pow(qpt[it][kp],2.0)/Math.pow(dpl [kp],5.333);
        Hwup[it][kp]=Hwdw[it][kp]+slopt[it][kp]*lp [kp];
       if(Hwup[it][kp]>=Hj[I0 [kp]])
		 {Hwup[it][kp]=Hj[I0 [kp]];
		  slopt[it][kp] = (Hwup[it][kp]-Hwdw[it][kp])/lp [kp];
		  vpt[it][kp]=Math.pow(rid[it]  [kp],0.6667)*Math.pow(slopt[it][kp],0.5)/slp [kp];
		  qpt [it][kp]=vpt[it][kp]*Ad0;
		 }}else
//--5
		 {hdc_min=(Hwdw[it][kp]-ZJdw [kp])/dpl [kp];
		  slopt[it][kp]=slop [kp];
		  if(hdc_min<0.0) hdc_min=0.0;
      ARD23=slp [kp]*qpt[it][kp]/Math.pow(slopt[it][kp],0.5);
//----------h/d---0.618aaa---------------------------
        hdc_max=1.0; fun_hd0=1000.0;
  	for(k1=1;k1<Nmax&&fun_hd0>0.0001;k1++)
//--6
		  {if(k1==1) {hdcc[0]=hdc_min; hdcc[1]=1.0;} else
		  {hdcc[0]=hdc_min+(hdc_max-hdc_min)*0.382;
   	   hdcc[1]=hdc_min+(hdc_max-hdc_min)*0.618;
		  }
         for(k00=0;k00<2;k00++)
		   {hda=1.0-2*hdcc[k00]; hdb=1.0-hdcc[k00]; hdc=Math.pow(Math.cos(hda),-1);
          hdd=hdcc[k00]*hdb;   hde=hda*Math.pow(hdd,0.5); hdf=Math.pow(dpl [kp],2.0);
          AD=hdf/4.0*hdc-hdf/2.00*hde;
		   RD=dpl [kp]/4.0 - dpl [kp]*hde/(2.0*hdc);
		   fun_hd[k00]=AD*Math.pow(RD,2.0/3.0) - ARD23;
		   }
     	   if(Math.abs(fun_hd[0])>Math.abs(fun_hd[1]))
		   {hdc_min=hdcc[0]; fun_hd0=Math.abs(fun_hd[0]);  hdcc0[it][kp]=hdcc[1];}
    	   else
		   {hdc_max=hdcc[1]; fun_hd0=Math.abs(fun_hd[1]);  hdcc0[it][kp]=hdcc[0];}
//
if(Iprt==1) /*outfile<<"it=  "<<it<<"  ik=  "<<ik<<"jk=  "<<jk<<"  kp=  "<<kp<<"  k1=  "<<k1<<"  hdcc[0]=  "<<hdcc[0]<<"  hdcc[1]=  "<<hdcc[1]
		 <<"  fun_hd[0]= "<<fun_hd[0]<<"  fun_hd[1]= "<<fun_hd[1]<<"  ARD23= "<<ARD23<<endl; */
		  //} //6-k1
	System.out.println("it=  "+it+"  ik=  "+ik+"jk=  "+jk+"  kp=  "+kp+"  k1=  "+k1+"  hdcc[0]=  "+hdcc[0]+"  hdcc[1]=  "+hdcc[1]
		 +"  fun_hd[0]= "+fun_hd[0]+"  fun_hd[1]= "+fun_hd[1]+"  ARD23= "+ARD23);
//
  rid [it][kp]=(0.25-0.5*hde/hdc)*dpl [kp];  
  Hwdw[it][kp]=ZJdw[kp]+hdcc0[it][kp]*dpl [kp];
  vpt[it][kp]=Math.pow(rid [it][kp],2.0/3.0)*Math.pow(slopt[it][kp],0.5)/slp [kp];
	Hwup[it][kp]=Hwdw[it][kp]+slopt[it][kp]*lp [kp];
//---------h/d---0.618aaa---------------------------
	 } //5--end_hdcc0[it][kp]
//
	/*outfile<<"   it= "<<it<<"   kp= "<<kp<<"   I0[kp]= "<<I0 [kp]<<"  Hwdm= "<<Hwdw[it][kp]<<"  Hj= "<<Hj[I0 [kp]]
		<<"  hdcc0= "<<hdcc0[it][kp]<<"  qpt= "<<qpt[it][kp]<<"  vpt= "<<vpt[it][kp]<<"  hdcc0= "<<hdcc0[it][kp]<<" fun_hd0= "<<fun_hd0<<endl; */
	System.out.println("   it= "+it+"   kp= "+kp+"   I0[kp]= "+I0 [kp]+"  Hwdm= "+Hwdw[it][kp]+"  Hj= "+Hj[I0 [kp]]
		+"  hdcc0= "+hdcc0[it][kp]+"  qpt= "+qpt[it][kp]+"  vpt= "+vpt[it][kp]+"  hdcc0= "+hdcc0[it][kp]+" fun_hd0= "+fun_hd0);
		 }//4 if(kp>=0) end
	}//3 ---jk end ---
	}//--2---ik end ---Hj[kp]
//---------------------------------  开始计算溢流节点  ---------------	
//
	 for(i=0;i<NP;i++)
	 {j=I0 [i]; 
	 if(Hwup [it][i]==Hj [j])
		{overflow [it][j]=overflow [it][j]-qpt [it][i];
	     for(ip=0;ip<NP;ip++)
			{k1=J0[ip];
	         if(k1==j) overflow [it][j]=overflow [it][j]+qpt [it][ip];
			}
       qjt [it][j]=Aj [j]*qit[it]*Acoef[j];
       overflow [it][j]=(overflow [it][j]+qjt [it][j])*dt*60.0;
		 Hw_over [it][j]=overflow [it][j]/Aj [j]/10000.0*1000.0;
		 }
	 }
//----qjt[it][j],overflow[it][j],totalflow[j],totalHw[j]
//------------------  计算溢流节点结束  ---------------	 
  //outfile<<endl;
  //outfile<<"    it   管段号  I0   J0   管径dpl   管段qp 水力半径R    充满度¨¨ 流速¨′(m/s)  上游水位  下游水位  上管底高  下管底高  上地面高"<<endl;
	System.out.println("\n    it   管段号  I0   J0   管径dpl   管段qp 水力半径R    充满度¨¨ 流速¨′(m/s)  上游水位  下游水位  上管底高  下管底高  上地面高");
	 for(i=0;i<NP;i++){
		//s.Format("%6d%6d%6d%5d%10.2lf%10.3lf%10.3lf%10.3lf%10.3lf%10.3lf%10.3lf%10.3lf%10.3lf%10.3lf\n",
		   //it,i,I0 [i],J0 [i],dpl [i],qpt [it][i],rid [it][i],hdcc0 [it][i],vpt [it][i],Hwup [it][i],Hwdw [it][i],ZJup [i],ZJdw [i],Hj [I0[i]]);
	  // outfile<<s;
		 System.out.println(it+i+I0 [i]+J0 [i]+dpl [i]+qpt [it][i]+rid [it][i]+hdcc0 [it][i]
				 +vpt [it][i]+Hwup [it][i]+Hwdw [it][i]+ZJup [i]+ZJdw [i]+Hj [I0[i]]);
  }
  //outfile<<endl;
	 System.out.println();
//}//1-- it end  ---
//------- 计算节点积水量和积水深度¨¨(m)----
	for(j=0;j<NN;j++) {totalflow [j]=0.0; totalHw [j]=0.0;}
  for(j=0;j<NN;j++)
	{for(it=0;it<NT;it++)
	{totalflow [j]=totalflow [j]+overflow [it][j];
	 totalHw [j]=totalHw [j]+Hw_over [it][j];
	} 
	}
//-----屏幕输出管段水力计算结束------
  //cout<<"    it   管段号  I0   J0   管径dpl   管段qp 水力半径R    充满度 流速(m/s)  上游水位  下游水位  上管底高  下管底高  上地面高”<<endl;
  System.out.println("    it   管段号  I0   J0   管径dpl   管段qp 水力半径R    充满度 流速(m/s)  上游水位  下游水位  上管底高  下管底高  上地面高");
//  
//  outfile<<endl; 
  for(it=0;it<NT;it++)
  {//outfile<<"  it= "<<it<<endl;
//  outfile<<endl;
//  outfile<<"    it   管段号  I0   J0   管径dpl   管段qp 水力半径R    充满度¨¨ 流速(m/s)  上游水位  下游水位    上管底高程 下管底高程<<endl;
	System.out.println("  it= ");
	System.out.println();
	System.out.println("    it   管段号  I0   J0   管径dpl   管段qp 水力半径R    充满度¨¨ 流速(m/s)  上游水位  下游水位    上管底高程 下管底高程");
	  for(i=0;i<NP;i++)
  {//s.Format("%6d%6d%6d%5d%10.2lf%10.3lf%10.3lf%10.3lf%10.3lf%10.3lf%10.3lf%10.3lf%10.3lf%10.3lf\n",
		  // it,i,I0 [i],J0 [i],dpl [i],qpt [it][i],rid [it][i],hdcc0 [it][i],vpt [it][i],Hwup [it][i],Hwdw [it][i],ZJup [i],ZJdw [i],Hj [I0[i]]);
//	   outfile<<s;
	  // cout<<s;
		System.out.println(it+i+I0 [i]+J0 [i]+dpl [i]+qpt [it][i]+rid [it][i]+hdcc0 [it][i]+vpt [it][i]+Hwup [it][i]+Hwdw [it][i]+ZJup [i]+ZJdw [i]+Hj [I0[i]]);
  }
//  outfile<<endl;
	}
//  outfile<<endl;
//-------------  节点溢流计算结果输出  ---------------	
  //outfile<<" ======== 时段节点积水量(m3) ========"<<endl;
  //outfile<<"  i=    ";
  System.out.println(" ======== 时段节点积水量(m3) ========");
  System.out.print("  i=    ");
	for(i=0;i<NN;i++) {if(i<10) {
		//outfile<<"  "<<i<<"   ";
		System.out.print("  "+i+"   ");
		} else  {
			//outfile<<" "<<i<<"   ";
			System.out.print(" "+i+"   ");
			}} 
	//outfile<<endl;
		System.out.println();
  //outfile<<"it="<<endl;
		System.out.println("it=");
	for(it=0;it<NT;it++)
  {if(it<10){
	  //outfile<<" "<<it<<"   ";
	  System.out.print(" "+it+"   ");
	  } else {
		  //outfile<<it<<"   ";
		  System.out.print(it + "   ");
	  }
  
//	 
	 for(i=0;i<NN;i++)
	 {if(overflow [it][i]<=0.0) {
		// outfile<<"      ";
		 System.out.print("      ");
		 } else
	 
	 { //s.Format("%6.1f",overflow [it][i]);
	   //outfile<<s;
			 System.out.print(overflow [it][i]);
	 }}
   //outfile<<endl;
	 System.out.println();
	}
//
  //outfile<<" ======== 时段节点积水深度输出(mm) ========"<<endl;
 // outfile<<"  i=    ";
	System.out.print(" ======== 时段节点积水深度输出(mm) ========");
	System.out.print("  i=    ");
	for(i=0;i<NN;i++) {if(i<10) {
		//outfile<<"  "<<i<<"   ";
		System.out.print("  "+i+"   ");
		} else  {
			//outfile<<" "<<i<<"   ";
			System.out.print(" "+i+"   ");
			}}  //outfile<<endl;
		System.out.println();	
 // outfile<<"it="<<endl;
		System.out.println("it=");
	for(it=0;it<NT;it++)
  {if(it<10){
	  //outfile<<" "<<it<<"   ";
	  System.out.print(" "+it+"   ");
	  } else {
		  //outfile<<it<<"   ";
		  System.out.print("   ");
	  }
  
//	 
	 for(i=0;i<NN;i++)
	 {if(overflow [it][i]<=0.0) {
		 //outfile<<"      ";
		 System.out.print("      ");
		 } else
	 
	 { //s.Format("%6.1f",Hw_over [it][i]);
	   //outfile<<s;
			 System.out.print(Hw_over [it][i]);
	 }}
   //outfile<<endl;
	 System.out.println();
	}
}//---000结束
//
	// outfile.close();
}}}
