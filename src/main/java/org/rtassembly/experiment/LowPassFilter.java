package org.rtassembly.experiment;

import java.io.BufferedReader;
import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;

import javax.sound.sampled.AudioFileFormat;
import javax.sound.sampled.AudioFormat;
import javax.sound.sampled.AudioInputStream;
import javax.sound.sampled.AudioSystem;
import javax.sound.sampled.DataLine;
import javax.sound.sampled.LineUnavailableException;
import javax.sound.sampled.TargetDataLine;

import org.jtransforms.fft.DoubleFFT_1D;
/*
 * This example records audio and saves the lowpass filtered data to a .wav file. 
 * There’s also an example of bandpass filtering. Filtering is done with FFT.
 * Note: this simple brick-wall filtering method is not applicable in some cases 
 * because it introduces ringing in the resulting filtered sound. 
 */
public class LowPassFilter implements Runnable {
    private final static int SAMPLERATE = 44100;
    private final static int BUFFERSIZE = SAMPLERATE * 2;

    private TargetDataLine tdl;
    private DoubleFFT_1D fft;

    LowPassFilter(TargetDataLine tdl) {
        this.tdl = tdl;
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

    private void lowPassFilter(double[] audioData, final int storedSamples) {
        // we need to initialize a buffer where we store our samples as complex numbers. first value is the real part, second is the imaginary.
        double[] fftData = new double[audioData.length * 2];
        for (int i = 0; i < storedSamples; i++) {
            // copying audio data to the fft data buffer, imaginary part is 0
            fftData[2 * i] = audioData[i];
            fftData[2 * i + 1] = 0;
        }

        // calculating the fft of the data, so we will have spectral power of each frequency component
        fft.complexForward(fftData);

        // zeroing out all components above this frequency
        int cutFreq = 2000;

        for (int i = 0; i < fftData.length; i += 2) {
            // lowpass
            if (i > ((cutFreq * (fftData.length/2)) / SAMPLERATE)*2)
                fftData[i] = fftData[i + 1] = 0;
        }
    
        // bandpass
        //if (i < (((cutFreq - 1000) * (fftData.length/2)) / SAMPLERATE)*2 || i > (((cutFreq + 1000) * (fftData.length/2)) / SAMPLERATE)*2)
            //fftData[i] = fftData[i + 1] = 0;

        // built-in scaling hangs the thread, so we don't use it
        fft.complexInverse(fftData, false);

        for (int i = 0; i < storedSamples; i++) {
            audioData[i] = fftData[2 * i] / (SAMPLERATE/4.0); // applying scaling
        }
    }

    @Override
    public void run() {
        byte[] abBuffer = new byte[tdl.getBufferSize()];
        double[] abBufferDouble = new double[abBuffer.length / 2];
        ByteArrayOutputStream baos = new ByteArrayOutputStream(); // this will store sound data

        fft = new DoubleFFT_1D(abBufferDouble.length);

        tdl.start();

        try {
            while (!Thread.interrupted()) {
                // waiting for the buffer to get filled
                while (tdl.available() < tdl.getBufferSize() * 0.5)
                    Thread.sleep(0, 1); // without this, the audio will be choppy

                int bytesRead = tdl.read(abBuffer, 0, tdl.available());

                // converting frames stored as bytes to double values
                int samplesRead = bytesRead / tdl.getFormat().getFrameSize();
                for (int i = 0; i < samplesRead; i++)
                    abBufferDouble[i] = ((abBuffer[i * 2] & 0xFF) | (abBuffer[i * 2 + 1] << 8)) / 32768.0;

                lowPassFilter(abBufferDouble, samplesRead);
                baos.write(getBytesFromDoubles(abBufferDouble, samplesRead), 0, samplesRead * 2);
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

