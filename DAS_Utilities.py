import shutil

def copyLargeFile(src, dest, buffer_size=1024*1024*1024):
    try:
        with open(src, 'rb') as fsrc:
            with open(dest, 'wb') as fdest:
                shutil.copyfileobj(fsrc, fdest, buffer_size)
    except:
        shutil.copyfile(src, dest)