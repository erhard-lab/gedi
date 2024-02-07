package executables;

import java.util.Locale;

import gedi.util.datastructure.collections.doublecollections.DoubleArrayList;
import gedi.util.math.stat.RandomNumbers;

public class SlamSimu {

	
	public static void main(String[] args) {
		
		
		double k_on=1.0/2;//(1.0/3); // inactive for 2h on average
		double k_off=1.0/(0.1); // active for 6 min on average
		double delta = 1; // half life of 1h
		double sigma = 10*k_off; // on average 10 mRNAs per burst
		
		double time = 0;
		int on = 0;
		DoubleArrayList mRNA = new DoubleArrayList();
		
		RandomNumbers rnd = new RandomNumbers();
		
		double nextStateChange = time+rnd.getExponential(on==1?k_off:k_on);
		double oldout = 100;
		for (int steps = 0; time<1000000; steps++) {
			if (time>1000&& time-1000>oldout) {
				int n = 0;
				int o = 0;
				for (int i=0; i<mRNA.size(); i++)
					if (mRNA.getDouble(i)<time-2)o++;
					else n++;
				System.out.printf(Locale.US,"%.3f\t%d\t%d\n",time,o,n);
				oldout = time;
			}
			if (on==1) {
				double syn = rnd.getExponential(sigma);
				double dec = rnd.getExponential(delta*mRNA.size());
				time=time+Math.min(syn, dec);
				if (time>nextStateChange) {
					time=nextStateChange;
					on = 1-on;
					nextStateChange = time+rnd.getExponential(on==1?k_off:k_on);
					continue;
				}
				
				if (syn<dec) 
					mRNA.add(time);
				else
					mRNA.removeEntry(rnd.getUnif(0, mRNA.size()));
			}
			else {
				double dec = rnd.getExponential(delta*mRNA.size());
				time=time+dec;
				if (time>nextStateChange) {
					time=nextStateChange;
					on = 1-on;
					nextStateChange = time+rnd.getExponential(on==1?k_off:k_on);
					continue;
				}
				mRNA.removeEntry(rnd.getUnif(0, mRNA.size()));
			}
			
		}
		
	}
	
	
}
