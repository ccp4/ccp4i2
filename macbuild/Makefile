#Make these sensible variables yourself, just examples.
TOP=$(shell (pwd))
MACEXENAME=ccp4i2
MACAPPNAME=ccp4i2.app
PROGNAME=ccp4i2
VERSION=15.3.13
OSNAME=macosx

$(MACEXENAME)-$(VERSION).dmg:
	 @(VOLNAME=$(MACEXENAME) ; \
	 mkdir -p MacOSX/BUILD ; \
	 DIRNAME=$(PROGNAME); \
	 TOP=$$PWD ; \
	 cd MacOSX/ ; \
	 sed -e "s/VERSIONNUMBER/$(VERSION)/" -e "s/CCP4MG/$(MACEXENAME)/g" Info.plist-template > Info.plist ; \
	 sed -e "s/VERSIONNUMBER/$(VERSION)/" -e "s/CCP4MG/$(MACEXENAME)/g" Resources/Description.plist-template > Resources/Description.plist ; \
	 cd BUILD ; \
	 rm -rf $$TOP/MacOSX/BUILD/$(MACAPPNAME) ; \
	 mkdir -p $(MACAPPNAME)/Contents ; \
	 cd $(MACAPPNAME)/Contents ; \
	 rm -fr ccp4mg ; \
	 rm -fr $(MACEXENAME) ; \
	 rm -fr $$DIRNAME ; \
	 tar zxf $$TOP/$(PROGNAME)-$(VERSION)-$(OSNAME).tar.gz ; \
	 mkdir -p Resources ; \
	 cp $$TOP/MacOSX/Resources/Description.plist Resources ; \
	 cp $$TOP/MacOSX/Resources/ccp4i2.icns Resources ; \
	 cp $$TOP/MacOSX/Info.plist . ; \
	 mkdir -p MacOS ; \
	 cd MacOS ; \
	 ln -s ../$$DIRNAME/bin/browser $(MACEXENAME) ; \
	 cd $$TOP/MacOSX/BUILD/ ; \
	 CCP4MGPKG=$$TOP/MacOSX/BUILD/$(MACAPPNAME) ; \
	 echo "Checking size of $$CCP4MGPKG" ; \
	 SIZE=`du -sk $$CCP4MGPKG | awk '{print $$1}'` ; \
	 SIZE=`expr  $$SIZE / 1000 + 200` ; \
	 cd $$TOP ; \
	 echo "Creating disk image $@" ; \
	 hdiutil create "$@" -megabytes "$$SIZE" -ov -type UDIF -quiet -fs HFS+ -volname "$$VOLNAME" ; \
	 echo "Attaching disk image $@" ; \
	 DEVICE=`hdiutil attach "$@" | awk '/Apple_HFS/{print $$1}'` ; \
	 DEVICE_AWK=`echo $$DEVICE | sed -e 's=/=\\\\/=g'` ; \
	 MOUNTPOINT=`df | awk '/'$$DEVICE_AWK'/{i=match($$0,"%");r=substr($$0,i);i=match(r,"/");r=substr(r,i);print r}'` ; \
	 echo "Copying $$CCP4MGPKG to $$MOUNTPOINT" ; \
	 cp -Rp "$$CCP4MGPKG" "$$MOUNTPOINT" ; \
	 ln -s /Applications "$$MOUNTPOINT/Applications" ; \
	 cp "$$TOP/MacOSX/Resources/simple_background.png" "$$MOUNTPOINT/background.png" ; \
	 /Developer/Tools/SetFile -a V "$$MOUNTPOINT/background.png" ; \
	 cp "$$TOP/MacOSX/Resources/ccp4i2.icns" "$$MOUNTPOINT/.VolumeIcon.icns" ; \
	 cp "$$TOP/MacOSX/QtMG-DS_Store" "$$MOUNTPOINT/.DS_Store" ; \
	 /Developer/Tools/SetFile -a C "$$MOUNTPOINT" ; \
	 echo "Ejecting $$MOUNTPOINT" ; \
	 hdiutil eject "$$DEVICE" -quiet ; \
	 hdiutil convert  "$@" -format UDZO -o "__$@" ; \
	 mv "__$@" $(MACEXENAME)-$(VERSION)-nolicence.dmg ; \
	 hdiutil convert  "$@" -format UDRO -o "_$@" ; \
	 hdiutil unflatten "_$@" ; \
	 /Developer/Tools/Rez /Developer/Headers/FlatCarbon/Carbon.r "$$TOP/MacOSX/SLA.r" -a -o "_$@" ; \
	 hdiutil flatten "_$@" ; \
	 hdiutil convert  "_$@" -format UDZO -o "__$@" ; \
	 mv "__$@" "$@" ; \
	 rm -f "_$@" ; \
	 echo "$@ is the final installation disk image" ; \
	 )

