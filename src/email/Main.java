package email;

public class Main {
	public static void main(String[] args) {
		String host = "smtp.163.com";            			// 设置邮件服务器
		String usrName = "m18508130880@163.com";	// 登录账号,一般都是和邮箱名一样吧
		String seneder = "m18508130880@163.com";  			// 发件人账号
		String pwd = "chenjian1234";        			// 发件人密码或授权码
		
		
		String recv = "903063674@qq.com";     			// 接收人账号
		String subject = "测试";            			// 主题
		String msg = "<h1>测试</h1>";     // 发送内容
		
		EmailUtil mailUtil = new EmailUtil();
		//初始化
		mailUtil.init(host, usrName, seneder, pwd, recv, subject, msg);
		//设置抄送人
		mailUtil.setCarbonCopy("3316281267@qq.com");
		mailUtil.setCarbonCopy("792750859@qq.com");
		mailUtil.send();
	}
}
