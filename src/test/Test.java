package test;

import java.io.File;

public class Test
{
	public static void main(String [] args){   
          
		String path = "C:/Users/Administrator/Desktop/雨水降雨模拟/数据表格";
        File file = new File(path);   
         
        File[] array = file.listFiles();   
          
        for(int i=0;i<array.length;i++){   
            if(array[i].isFile()){   
                // only take file name   
                System.out.println("^^^^^" + array[i].getName());   
                // take file path and name   
                //System.out.println("#####" + array[i]);   
                // take file path and name   
                //System.out.println("*****" + array[i].getPath());   
            }else if(array[i].isDirectory()){   
                //getFile(array[i].getPath()); 
            	System.out.println("失败！");
            }   
        }  
	}
}
