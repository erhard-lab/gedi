package gedi.util.r;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.lang.ProcessBuilder.Redirect;
import java.net.ConnectException;
import java.net.ServerSocket;
import java.net.Socket;
import java.net.SocketException;
import java.net.UnknownHostException;
import java.util.concurrent.atomic.AtomicInteger;

import cern.colt.Arrays;
import gedi.util.StringUtils;
import gedi.util.functions.EI;
import gedi.util.functions.ParallelizedState;

public class RProcess implements AutoCloseable {

	private Thread thread;
	private Process process;
	private int port = -1;
	private BufferedReader rReader;
	private BufferedWriter rWriter;
	private Socket socket;
	
	// Thread-local instance storage
    private static final ThreadLocal<RProcess> threadLocalInstance = ThreadLocal.withInitial(() -> {
        try {
            return new RProcess(false);
        } catch (IOException e) {
            throw new RuntimeException("Failed to start R process", e);
        }
    });
    
    // Public accessor
    public static RProcess getForCurrentThread() {
        return threadLocalInstance.get();
    }
	
	public RProcess(boolean interactive) throws IOException {
		ProcessBuilder b = new ProcessBuilder();
		b.redirectError(Redirect.INHERIT);
		if (interactive) {
			b.redirectOutput(Redirect.INHERIT);
			b.command("R", "--interactive","--no-save", "--no-restore");
		} else {
	        b.command("R", "--slave", "--vanilla");
		}
		
		process = b.start();
		
		thread = new Thread(()->{
			try {
				process.waitFor();
			} catch (Exception e) {
			} 
		});
		thread.setDaemon(true);
		thread.setName("R-Process");
		thread.start();
		
		rReader = new BufferedReader(new InputStreamReader(process.getInputStream()));
		rWriter = new BufferedWriter(new OutputStreamWriter(process.getOutputStream()));
	   
	}
	
	public BufferedWriter getStdin() {
		return rWriter;
	}
	
	public BufferedReader getStdout() {
		return rReader;
	}
	
	public void read() throws IOException {
		read(System.in);
	}
	
	public void read(InputStream is) throws IOException {
		int key;
		InputStreamReader in = new InputStreamReader(is);
		while (process.isAlive()) {
			key = in.read();
			if (key==10)
				return;
			if (process.isAlive()) {
				getStdin().write(key);
				getStdin().flush();
			}
		}
	}
    @Override
    public void close() {
        try {
            synchronized (rWriter) {
                rWriter.write("quit(save='no')\n");
                rWriter.flush();
            }
            rWriter.close();
            rReader.close();
            socket.close();
            process.waitFor();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

	
	public static final int start_port = 12844;
	public static final AtomicInteger counter = new AtomicInteger();
	public boolean isRunning() {
		return process.isAlive();
	}
	
	private synchronized void connectIfNotConnected() throws IOException {
		if (port==-1) {
			ServerSocket probe = new ServerSocket(0);
			port = probe.getLocalPort();
			probe.close();
			getStdin().append("con<-socketConnection(port = "+port+", blocking = TRUE, server = TRUE, open = \"rb\")\n").flush();
			while (socket==null) {
				try {
					socket = new Socket("localhost",port);
				} catch (ConnectException e) {
				}
			}
		}
	}
	
	private synchronized void acceptConnection() throws IOException {
		if (port==-1) throw new RuntimeException("Not connected!");
		getStdin().append(".Internal(loadFromConn2(con,environment(),F))\n").flush();
	}

	public RDataWriter startSetting() throws UnknownHostException, IOException {
		 connectIfNotConnected(); 
		 acceptConnection();
		
		RDataWriter re = new RDataWriter(socket.getOutputStream());
//		re.setFinishCallback(()->{
//			try {
//				eval("close(con)");
//			} catch (IOException e) {
//				throw new RuntimeException(e);
//			}
//		});
		re.writeHeader();
		return re;
	}
	
	public static String getListenCommand(int port) {
//		return ".Internal(loadFromConn2(socketConnection(port="+port+",blocking=T,server=T,open=\"r+b\"),environment(),F))\n";
		return "con <- socketConnection(port = "+port+", blocking = TRUE, server = TRUE, open = \"rb\")\n"
				+ ".Internal(loadFromConn2(con,environment(),F))\n"
				+ "close(con)\n";
	}
	
	public String eval(String cmd) throws IOException {
		String wrappedCmd = String.format("cat(paste(%s, collapse=','), '\\n')\n", cmd);
        getStdin().append(wrappedCmd).flush();
        // Read single line result
        String line;
        while ((line = rReader.readLine()) != null) {
            if (!line.trim().isEmpty()) {
                return line;
            }
        }
        throw new IOException("No output from R process");
	}
	
	public void library(String library) throws IOException {
		eval("library("+library+")");
	}
	

	public double[] callNumericFunction(String cmd) throws IOException {
       String re = eval(cmd);
       return StringUtils.parseDouble(StringUtils.split(re,','));       
	}
	
	public static void main(String[] args) throws IOException {
		
		RProcess proc = RProcess.getForCurrentThread();
		proc.startSetting().write("a", new double[] {0,0}).finish();
		double[] a = proc.callNumericFunction("a+1");
		System.out.println(Arrays.toString(a));

		RDataWriter setter = proc.startSetting();
		
		setter.write("ll0", new double[] {-0.984425711460142});
		setter.write("ll1", new double[] {-0.07074455952965866});
		setter.write("c", new int[] {1});
		setter.finish();
		
		proc.library("grandR");
		double[] re = proc.callNumericFunction("grandR:::fit.ntr.betamix(ll0,ll1,c)");
		
		System.out.println(Arrays.toString(re));
			
		
		
		EI.seq(0, 10).parallelized(5, 1, ei->ei.map(ind -> {
			try {
				RProcess rproc = RProcess.getForCurrentThread();
				rproc.startSetting().write("a", new double[] {ind,ind*2}).finish();
				double[] aa = rproc.callNumericFunction("a+1");
				return aa[0];
			} catch (Exception e) {
				throw new RuntimeException(e);
			}
		})).sort().print();
	}

	
}
