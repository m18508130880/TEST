package test;

public class M
{

	public static void main(String[] args)
	{

		int[][] a = { { 1, 2, 3 }, { 4, 5, 6 } };
		// 定义转换矩阵数组
		int[][] b = new int[a[0].length][a.length];
		// 给转换矩阵数组赋值
		for (int i = 0; i < a.length; i++)
		{
			for (int j = 0; j < a[0].length; j++)
			{
				b[j][i] = a[i][j];
			}
		}

		// 输出转换后的矩阵数组
		for (int i = 0; i < b.length; i++)
		{
			for (int j = 0; j < b[0].length; j++)
			{
				System.out.print(b[i][j] + " ");
			}
			System.out.println();
		}

	}

}