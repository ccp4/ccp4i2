#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <CoreFoundation/CoreFoundation.h>

int main (int argc, char **argvin) {
    CFBundleRef mainBundle;
    CFURLRef bundleURL;
    mainBundle = CFBundleGetMainBundle();
    bundleURL = CFBundleCopyBundleURL(mainBundle);
    CFStringRef str = CFURLCopyFileSystemPath( bundleURL, kCFURLPOSIXPathStyle );
    CFRelease(bundleURL);
    char path[PATH_MAX];
    CFStringGetCString( str, path, FILENAME_MAX, kCFStringEncodingASCII );
    CFRelease(str);
    strncat(path,"/../bin/ccp4i2",14);
    printf("%s\n",path);
    char * const argv[] = {path, NULL };
    if(execv(path, argv) < 0) {
      perror("execv error");
    }
    return 0;
}
