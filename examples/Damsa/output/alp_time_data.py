import numpy as np
import matplotlib.pyplot as plt

# 데이터 파일 읽기
data = np.loadtxt('output/alp_time_data.txt')

# 히스토그램 그리기
plt.hist(data, bins=50, alpha=0.75, color='blue', edgecolor='black')

# 그래프 제목과 레이블 설정
plt.title('Distribution of Alp Time Data')
plt.xlabel('Value')
plt.ylabel('Frequency')

# 그래프 저장하기
plt.grid(True)
plt.savefig('output/alp_time_histogram.png')  # 그래프를 PNG 파일로 저장
plt.close()  # 그래프 창 닫기
