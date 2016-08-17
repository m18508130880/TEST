package email;

import java.util.ArrayList;

import org.apache.commons.mail.HtmlEmail;    
    
public class EmailUtil { 
	
	//初始化
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
	
	//设置抄送
	public void setCarbonCopy(String cc){
		if(cc != null){
			carbonCopy.add(cc);
		}
	}
	
	// 发送email  
    public void send() {  
        HtmlEmail email = new HtmlEmail();  
        try {  
        	// 字符编码集的设置  
        	email.setCharset("UTF-8");  
            // 这里是SMTP发送服务器的名字：163的如下："smtp.163.com"
        	// 要发送邮件的邮箱必须开通SMTP协议
            email.setHostName(host);  
            // 发送人的邮箱  
            email.setFrom(sender, "");  
            // 如果需要认证信息的话，设置认证：用户名-密码。分别为发件人在邮件服务器上的注册名称和密码  
            email.setAuthentication(username, password);  
            // 收件人的邮箱  
            email.addTo(receiver);
            //抄送人邮箱
            if(carbonCopy != null){
    			for(int i = 0;i < carbonCopy.size();i ++){
    				email.addCc(carbonCopy.get(i));
    			}
    		}
            // 要发送的邮件主题  
            email.setSubject(subject);  
            // 要发送的信息，由于使用了HtmlEmail，可以在邮件内容中使用HTML标签  
            email.setMsg(message);
            // 发送  
            email.send();  
            System.out.println("发送邮件成功！");
        } catch (Exception e) {  
        	System.out.println("发送邮件失败！");
        }  
    }
    private String host;     // 服务器地址  
    
    private String sender;   // 发件人的邮箱  
  
    private String receiver; // 收件人的邮箱
    
    private ArrayList<String> carbonCopy; // 抄送
    
	private String username; // 账号  
  
    private String password; // 密码  
  
    private String subject;  // 主题  
  
    private String message;  // 信息(支持HTML) 
    
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