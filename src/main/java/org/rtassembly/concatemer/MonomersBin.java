package org.rtassembly.concatemer;

import java.util.HashMap;
import java.util.List;

public class MonomersBin {
	private static int _ID_=0;
	
	public int id;
	short[] monomer;
	int count;
	HashMap<String, short[]> monomersList;
	
	public MonomersBin(){
		id=_ID_++;
	}
	public MonomersBin(short[] mnm){
		this();
		monomersList=new HashMap<>();
		monomer=mnm;
	}
	
	public int getLength(){
		return monomer==null?0:monomer.length;
	}
	
	public int getCount(){
		return count;
	}

	public short[] getMonomerSignal(){
		return monomer;
	}
	
	//chop the concatemeric signal at breakPoints into monomers and append to this bin
	public void addToBin(String rname, short[] signal, List<Integer> breaks, int lag){
		
		
		flushToFile();
	}
	//do this every time the bin reaches 100 monomers
	void flushToFile(){
		String fname="bin_"+id;
		
	}
}
