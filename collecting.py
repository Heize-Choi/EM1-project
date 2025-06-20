import cv2
import datetime
import csv
import time
import schedule
import logging
import sys

def measure_illuminance(writer):
    try:
        ret, frame = cap.read()
        if not ret:
            logging.warning("카메라 프레임 읽기 실패")
            return

        gray = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)
        illuminance = gray.mean()
        timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        writer.writerow([timestamp, f"{illuminance:.2f}"])
        logging.info(f"{timestamp} - Illuminance: {illuminance:.2f}")
    except Exception as e:
        logging.error(f"측정 중 예외: {e}")

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s [%(levelname)s] %(message)s",
                        handlers=[logging.StreamHandler(sys.stdout)])
    # 종료 시각 입력
    end_time_str = input("종료 시각을 입력하세요 (YYYY-MM-DD HH:MM:SS): ")
    try:
        end_time = datetime.datetime.strptime(end_time_str, "%Y-%m-%d %H:%M:%S")
    except ValueError:
        print("날짜 형식이 잘못되었습니다.")
        sys.exit(1)

    # CSV 파일 준비
    f = open('illuminance4.csv', 'w', newline='')
    writer = csv.writer(f)
    writer.writerow(['Timestamp', 'Illuminance'])

    # 카메라 열기
    cap = cv2.VideoCapture(0)
    if not cap.isOpened():
        print("카메라를 열 수 없습니다.")
        f.close()
        sys.exit(1)

    # 스케줄 등록 (10분마다)
    schedule.every(10).minutes.do(measure_illuminance, writer=writer)

    # 측정 루프
    try:
        while True:
            now = datetime.datetime.now()
            if now >= end_time:
                print("종료")
                break
            schedule.run_pending()
            time.sleep(1)
    finally:
        cap.release()
        f.close()
        print("프로그램 종료.")
