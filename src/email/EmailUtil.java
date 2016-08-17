package email;

import java.util.ArrayList;

import org.apache.commons.mail.HtmlEmail;    
    
public class EmailUtil { 
	
	//��ʼ��
	public void init(String host,String usrName,String sender,String pwd,String recv,String subject,String msg){
		this.host = host;
		this.username = usrName;
		this.sender = sender;
		this.password = pwd;
		
		this.receiver = recv;
		carbonCopy = new ArrayList();
		this.subject = subject;
		this.message = msg;
	}
	
	//���ó���
	public void setCarbonCopy(String cc){
		if(cc != null){
			carbonCopy.add(cc);
		}
	}
	
	// ����email  
    public void send() {  
        HtmlEmail email = new HtmlEmail();  
        try {  
        	// �ַ����뼯������  
        	email.setCharset("UTF-8");  
            // ������SMTP���ͷ����������֣�163�����£�"smtp.163.com"
        	// Ҫ�����ʼ���������뿪ͨSMTPЭ��
            email.setHostName(host);  
            // �����˵�����  
            email.setFrom(sender, "");  
            // �����Ҫ��֤��Ϣ�Ļ���������֤���û���-���롣�ֱ�Ϊ���������ʼ��������ϵ�ע�����ƺ�����  
            email.setAuthentication(username, password);  
            // �ռ��˵�����  
            email.addTo(receiver);
            //����������
            if(carbonCopy != null){
    			for(int i = 0;i < carbonCopy.size();i ++){
    				email.addCc(carbonCopy.get(i));
    			}
    		}
            // Ҫ���͵��ʼ�����  
            email.setSubject(subject);  
            // Ҫ���͵���Ϣ������ʹ����HtmlEmail���������ʼ�������ʹ��HTML��ǩ  
            email.setMsg(message);
            // ����  
            email.send();  
            System.out.println("�����ʼ��ɹ���");
        } catch (Exception e) {  
        	System.out.println("�����ʼ�ʧ�ܣ�");
        }  
    }
    private String host;     // ��������ַ  
    
    private String sender;   // �����˵�����  
  
    private String receiver; // �ռ��˵�����
    
    private ArrayList<String> carbonCopy; // ����
    
	private String username; // �˺�  
  
    private String password; // ����  
  
    private String subject;  // ����  
  
    private String message;  // ��Ϣ(֧��HTML) 
    
	public String getHost() {  
        return host;  
    }  
  
    public void setHost(String host) {  
        this.host = host;  
    }  
  
    public String getSender() {  
        return sender;  
    }  
  
    public void setSender(String sender) {  
        this.sender = sender;  
    }  
  
    public String getReceiver() {  
        return receiver;  
    }  
  
    public void setReceiver(String receiver) {  
        this.receiver = receiver;  
    }  
  
    public String getUsername() {  
        return username;  
    }  
  
    public void setUsername(String username) {  
        this.username = username;  
    }  
  
    public String getPassword() {  
        return password;  
    }  
  
    public void setPassword(String password) {  
        this.password = password;  
    }  
  
    public String getSubject() {  
        return subject;  
    }  
  
    public void setSubject(String subject) {  
        this.subject = subject;  
    }  
  
    public String getMessage() {  
        return message;  
    }  
  
    public void setMessage(String message) {  
        this.message = message;  
    }
  
}  