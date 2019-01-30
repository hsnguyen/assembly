/*	This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>. */

package org.rtassembly.experiment;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.io.UnsupportedEncodingException;

import javax.sound.sampled.AudioFileFormat;
import javax.sound.sampled.AudioFormat;
import javax.sound.sampled.AudioInputStream;
import javax.sound.sampled.AudioSystem;
import javax.sound.sampled.DataLine;
import javax.sound.sampled.LineUnavailableException;
import javax.sound.sampled.TargetDataLine;

import org.jtransforms.fft.DoubleFFT_1D;

/* 
 * Lowpass FIR with FFT convolution to eliminate ring effects
 * http://www.dspguide.com/ch18/2.htm
 */

public class LowPassFilter2 implements Runnable {
	private final static int SAMPLERATE = 48000;
	private final static int BUFFERSIZE = SAMPLERATE * 2;
	private final static int FFTRES = 512;

	private TargetDataLine tdl;
	private DoubleFFT_1D fft;
	private byte[] audioBytes;
	private ByteArrayOutputStream baos; // this will store sound data
	private double[] filter, filterFFT, audio, audioFFT, audioOverlap, audioFiltered;

	// see http://en.wikipedia.org/wiki/Sinc_function
	private double sinc(double x) {
		if (x == 0)
			return 1;
		return Math.sin(Math.PI*x)/(Math.PI*x);
	}
	
	LowPassFilter2(TargetDataLine tdl) {
		audio = new double[FFTRES/2];
		filter = new double[FFTRES/2+1];
		audioFFT = new double[FFTRES*2]; // *2 because FFT data has both real and imaginary parts
		filterFFT = new double[FFTRES*2];
		audioOverlap = new double[FFTRES/2];
		audioFiltered = new double[FFTRES/2];
		
		this.tdl = tdl;
		audioBytes = new byte[audio.length*2];
		baos = new ByteArrayOutputStream();
		fft = new DoubleFFT_1D(FFTRES);
		System.out.println("FFT resolution: " + FFTRES);
		System.out.println("Audio sample in buffer length: " + audio.length);
		System.out.println("Filter length: " + filter.length);

		// designing the windowed sinc filter
		int cutFreq = 1000;
		for (int i = 0; i < filter.length-1; i++) {
			// see http://en.wikipedia.org/wiki/Sinc_filter
			double sincFilter = (2*cutFreq/(double)SAMPLERATE)*sinc(2*cutFreq*(i-((filter.length-1)/2.0))/(double)SAMPLERATE);

			// uncomment for a bandpass filter
			//sincFilter = sincFilter-(2*(cutFreq+1000)/(double)SAMPLERATE)*sinc(2*(cutFreq+1000)*(i-((filter.length-1)/2.0))/(double)SAMPLERATE);

			// uncomment for a highpass filter
			/*if (i != (filter.length-1)/2)
				sincFilter *= -1;
			else
				sincFilter = 1-sincFilter;*/

			// applying a Hamming window, see http://stackoverflow.com/questions/5418951/what-is-the-hamming-window-for and http://en.wikipedia.org/wiki/Window_function
			//filter[i] = (0.53836 - (0.46164 * Math.cos((Math.PI*2) * (double)i / (double)(filter.length-1)))) * idealFilter;
			// applying a Hann window
			//filter[i] = 0.5*(1-Math.cos( (2*Math.PI*i)/(double)(filter.length-1) )) * sincFilter;
			// applying a Blackman window
			filter[i] = (0.42-0.5*Math.cos((2*Math.PI*i)/(double)(filter.length-1))+0.08*Math.cos((4*Math.PI*i)/(double)(filter.length-1))) * sincFilter;
		}

		//baos.write(getBytesFromDoubles(filter, filter.length), 0, filter.length * 2); writeWavFile(baos.toByteArray(), baos.size() / 2, "output.wav"); System.exit(0);

		// clearing the audio overlap buffer
		for (int i = 0; i < audioOverlap.length; i++)
			audioOverlap[i] = 0;
		// clearing the FFT buffer
		for (int i = 0; i < filterFFT.length; i++)
			filterFFT[i] = 0;
		// copying time domain filter data to the FFT buffer
		for (int i = 0; i < filter.length; i++)
			filterFFT[2 * i] = filter[i];

		/*try {
			BufferedWriter outputStream = new BufferedWriter(new OutputStreamWriter(new FileOutputStream("output.txt"), "UTF-8"));
			for (int i = 0; i < filterFFT.length; i++)
				outputStream.write(i + " "+filterFFT[i]+"\n");
			outputStream.close();
		} catch (Exception e) {
			e.printStackTrace();
		}*/

		fft.complexForward(filterFFT);
	}

	// converts float array to byte array
	private byte[] getBytesFromDoubles(final double[] audioData, final int storedSamples) {
		byte[] audioDataBytes = new byte[storedSamples * 2];

		for (int i = 0; i < storedSamples; i++) {
			// saturation
			audioData[i] = Math.min(1.0, Math.max(-1.0, audioData[i]));

			// scaling and conversion to integer
			int sample = (int) Math.round((audioData[i] + 1.0) * 32767.5) - 32768;

			byte high = (byte) ((sample >> 8) & 0xFF);
			byte low = (byte) (sample & 0xFF);
			audioDataBytes[i * 2] = low;
			audioDataBytes[i * 2 + 1] = high;
		}

		return audioDataBytes;
	}

	// saves the audio data given in audioDataBytes to a .wav file
	private void writeWavFile(final byte[] audioDataBytes, final int storedSamples, final String fileName) {
		AudioFormat audioFormat = new AudioFormat(AudioFormat.Encoding.PCM_SIGNED, SAMPLERATE, 16, 1, 2, SAMPLERATE, false);
		AudioInputStream audioInputStream = new AudioInputStream(new ByteArrayInputStream(audioDataBytes), audioFormat, storedSamples);

		try {
			FileOutputStream fileOutputStream = new FileOutputStream(fileName);
			AudioSystem.write(audioInputStream, AudioFileFormat.Type.WAVE, fileOutputStream);
			audioInputStream.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	private void lowPassFilter(double[] audioDataIn, double[] audioDataOut) {
		// zeroing out the audio FFT buffer
		for (int i = 0; i < audioFFT.length; i++)
			audioFFT[i] = 0;
		// copying audio data to the FFT data buffer
		for (int i = 0; i < audioDataIn.length; i++)
			audioFFT[2 * i] = audioDataIn[i];

		// calculating the fft of the data
		fft.complexForward(audioFFT);

		// pointwise multiplication of the filter and audio data in the frequency domain
		for (int i = 0; i < filterFFT.length; i += 2) {
			double temp = audioFFT[i] * filterFFT[i] - audioFFT[i+1]*filterFFT[i+1];
			audioFFT[i+1] = audioFFT[i] * filterFFT[i+1] + audioFFT[i+1] * filterFFT[i]; // imaginary part
			audioFFT[i] = temp; // real part
		}

		// built-in scaling hangs the thread, so we don't use it
		fft.complexInverse(audioFFT, false);

		// adding the first half of the audio FFT buffer to the overlap buffer
		for (int i = 0; i < audioDataOut.length; i++)
			audioDataOut[i] = (audioOverlap[i] + audioFFT[i * 2]) / 2000; // applying scaling

		// copying the second half of the audio FFT buffer to the audio overlap buffer
		for (int i = 0; i < audioOverlap.length; i++)
			audioOverlap[i] = audioFFT[audioFFT.length/2 + i * 2];
	}

	@Override
	public void run() {
		tdl.start();

		try {
			while (!Thread.interrupted()) {
				// waiting for the buffer to get filled
				while (tdl.available() < audioBytes.length)
					Thread.sleep(0, 1); // without this, the audio will be choppy

				int	bytesRead = tdl.read(audioBytes, 0, audioBytes.length);

				// converting frames stored as bytes to double values
				int samplesRead = bytesRead / 2;
				for (int i = 0; i < samplesRead; i++)
					audio[i] = ((audioBytes[i * 2] & 0xFF) | (audioBytes[i * 2 + 1] << 8)) / 32768.0;

				lowPassFilter(audio, audioFiltered);
				baos.write(getBytesFromDoubles(audioFiltered, audioFiltered.length), 0, audioFiltered.length * 2);
			}
		} catch (InterruptedException e) {
		}

		tdl.stop();
		tdl.close();

		writeWavFile(baos.toByteArray(), baos.size() / 2, "output.wav");
	}

	public static void main(String[] args) {
		AudioFormat audioFormat = new AudioFormat(AudioFormat.Encoding.PCM_SIGNED, SAMPLERATE, 16, 1, 2, SAMPLERATE, false);
		DataLine.Info info = new DataLine.Info(TargetDataLine.class, audioFormat, BUFFERSIZE);

		TargetDataLine targetDataLine = null;
		try {
			targetDataLine = (TargetDataLine) AudioSystem.getLine(info);
			targetDataLine.open(audioFormat, BUFFERSIZE);
			System.out.println("Buffer size: " + targetDataLine.getBufferSize());
		} catch (LineUnavailableException e1) {
			e1.printStackTrace();
		}

		// creating the recorder thread from this class' instance
		LowPassFilter lowPassFilter = new LowPassFilter(targetDataLine);
		Thread lowPassFilterThread = new Thread(lowPassFilter);

		// we use this to read line from the standard input
		BufferedReader br = new BufferedReader(new InputStreamReader(System.in));
		lowPassFilterThread.setPriority(Thread.MAX_PRIORITY);
		lowPassFilterThread.start();

		System.out.println("Recording... press ENTER to stop recording!");
		try {
			br.readLine();
		} catch (IOException e) {
			e.printStackTrace();
		}

		lowPassFilterThread.interrupt();

		try {
			// waiting for the recorder thread to stop
			lowPassFilterThread.join();
		} catch (InterruptedException e) {
			e.printStackTrace();
		}
		System.out.println("Recording stopped.");
	}
}