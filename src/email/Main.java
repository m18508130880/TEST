package email;

public class Main {
	public static void main(String[] args) {
		String host = "smtp.163.com";            			// �����ʼ�������
		String usrName = "m18508130880@163.com";	// ��¼�˺�,һ�㶼�Ǻ�������һ����
		String seneder = "m18508130880@163.com";  			// �������˺�
		String pwd = "chenjian1234";        			// �������������Ȩ��
		
		
		String recv = "903063674@qq.com";     			// �������˺�
		String subject = "����";            			// ����
		String msg = "<h1>����</h1>";     // ��������
		
		EmailUtil mailUtil = new EmailUtil();
		//��ʼ��
		mailUtil.init(host, usrName, seneder, pwd, recv, subject, msg);
		//���ó�����
		mailUtil.setCarbonCopy("3316281267@qq.com");
		mailUtil.setCarbonCopy("792750859@qq.com");
		mailUtil.send();
	}
}
