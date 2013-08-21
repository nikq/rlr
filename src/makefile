# Makefile

# プログラム名とオブジェクトファイル名
program = rlr
objs    = scene.o vectormath.o modelholder.o mtseq.o

# 定義済マクロの再定義
CC     = g++ #gcc
CXX    = g++
CFLAGS = -Wall -O2 -msse -funroll-loops

# サフィックスルール適用対象の拡張子の定義
.SUFFIXES: .cxx .o

# プライマリターゲット
$(program): $(objs)
	$(CC) -o $(program) $^

# サフィックスルール
.cxx.o:
	$(CXX) $(CFLAGS) -c $<

# ファイル削除用ターゲット
.PHONY: clean
clean:
	$(RM) $(program) $(objs)


scene.o       : scene.cxx scene.h bvh.h vectormath.h modelholder.h kdtree.h
vectormath.o  : vectormath.cxx vectormath.h
modelholder.o : modelholder.cxx modelholder.h
mtseq.o       : mtseq.cxx mtseq.h
