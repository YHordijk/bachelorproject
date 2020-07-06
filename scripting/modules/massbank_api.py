import urllib.request
import numpy as np


def get_peaks(ID):
	url = f"http://massbank.jp/RecordDisplay?id={ID}"
	with urllib.request.urlopen(url) as site:
		lines = site.readlines()
		start_peaks = False
		peaks = []
		for line in lines:
			if start_peaks:
				if b'//\n' in line:
					start_peaks = False
					break
				peaks.append(line)
			if b'PK$PEAK' in line:
				start_peaks = True
				

	if len(peaks) == 0: return

	new_peaks = []
	for peak in peaks:
		peak = peak.decode("utf-8")
		peak = peak.strip('\n')
		peak = peak.replace('&nbsp', ' ')
		peak = peak.strip(' ')
		peak = peak.split(' ')

		new_peaks.append((float(peak[0]), float(peak[1])))

	#normalize peaks
	new_peaks = [(p[0], p[1]/np.asarray(new_peaks)[:,1].max()) for p in new_peaks]
	return new_peaks


def save_peaks(file, peaks):
	with open(file, 'w+') as f:
		for peak in peaks:
			f.write(str(peak[0]) + ', ' + str(peak[1]) + '\n')


def load_peaks(file):
	peaks = []
	with open(file, 'r') as f:
		for line in f.readlines():
			peak = line.split(',')
			peak = str(peak[0]).strip(' '), str(peak[1]).strip(' ')
			peak = float(peak[0]), float(peak[1])
			peaks.append(peak)

	return peaks


def get_mol(ID):
	url = f"http://massbank.jp/RecordDisplay?id={ID}"
	with urllib.request.urlopen(url) as site:
		lines = site.readlines()
		for line in lines:
			if b'CH$NAME' in line:
				line = line.decode('utf-8')
				line = line.strip('\n')
				line = line.replace('</b>', '')
				line = line.split(':')[1]
				line = line.strip('<br>')
				line = line.strip(' ')
				return line


if __name__ == '__main__':
	print(get_peaks(ID="TY000256"))
