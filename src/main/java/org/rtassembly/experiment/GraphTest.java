package org.rtassembly.experiment;

import java.io.BufferedReader;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import japsa.seq.SequenceReader;

public class GraphTest {
    public static void main(String[] args) throws Exception {
        // creating the process 
        ProcessBuilder pb = new ProcessBuilder("bwa").redirectErrorStream(true); 
          
        // map view of this process builder's environment 
        Map<String, String> envMap = pb.environment(); 
          
        // checking map view of environment 
        for(Map.Entry<String, String> entry : envMap.entrySet()){
             // checking key and value separately 
            System.out.println("Key = " + entry.getKey() +  
               ", Value = " + entry.getValue()); 
        } 
        System.out.println("System:\n" + System.getenv("PATH"));
        
		Process process =  pb.start();
		BufferedReader bf = SequenceReader.openInputStream(process.getInputStream());


		String line;
		String version = "";
		Pattern versionPattern = Pattern.compile("^Version:\\s(\\d+\\.\\d+\\.\\d+).*");
		Matcher matcher=versionPattern.matcher("");
		
		while ((line = bf.readLine())!=null){				
			matcher.reset(line);
			if (matcher.find()){
			    version = matcher.group(1);
			    System.out.println(line + " => version=" + version);
			    break;//while
			}
			
							
		}	
		bf.close();

    	
    }
}
