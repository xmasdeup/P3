import numpy as np
import matplotlib.pyplot as plt
import os
import glob
import wave

file = 'pitch_db/train/rl006.wav'


with  wave.open(file,'rb') as wav:
    data = wav.readframes(wav.getnframes())
    rate = wav.getframerate()

length = np.round(rate * 0.03,0)

data_array = np.frombuffer(data, dtype=np.int16)

index = np.argmax(np.abs(data_array))

shift = np.int16(index + length)

frame = data_array[index:shift]

correl = []

for i in enumerate(frame):

    correl.append( np.int32(frame[i[0]]) * np.int32(frame[i[0]]))

correl = np.array(correl)

x = np.linspace(0,1,len(frame))

plt.subplot(2,1,1)
plt.plot(x,frame, label= 'Mostra de veu')

plt.subplot(2,1,2)
plt.plot(x,correl, label= 'Correlació de la mostra')

plt.tight_layout()

plt.show()