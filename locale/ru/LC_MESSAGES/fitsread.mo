��    U      �  q   l      0     1  '   F  /   n     �     �     �     �  '   �          (     B      Q     r     �     �     �     �     �     �     	     ,	  
   J	  &   U	  !   |	      �	  $   �	  -   �	  !   
  ,   4
     a
     z
     �
     �
     �
     �
     �
  .        0  -   F  &   t     �  "   �     �  5   �  +        A     X  U   p  %   �     �               2  #   F     j      �     �     �  &   �  (   �  %   %  '   K     s     �     �  
   �     �     �  F   �             *   .  9   Y     �  �   �  3   F     z     �     �     �     �  <   �  #   ;     _  A  n     �  (   �  ;   �      -     N     m     �  C   �     �            .   +  !   Z     |     �     �     �     �  $        ,     K     e  +   q     �     �  1   �  ?   	  +   I  )   u  $   �     �     �               !     5  6   O     �  5   �  6   �     	  5   %     [  <   u  &   �     �     �  I     )   [  !   �  !   �  #   �     �       &   !     H     a     ~  ?   �  0   �  .     )   4     ^     u     �     �     �     �  ^   �     2     @  ;   O  C   �      �  �   �  D   �     �  '   �  '     #   9  %   ]  N   �  7   �     
         '   K          6           <       :          #      I   "   +   P   	       L                      )          7          H   2       ?   /   ;   $       8       .   (   @   J         >   R      U   N   E       A       O   5       =      3                     0   F   B   1   &   -          C                           !      G   *      Q          D          %           S   9   4   ,       M      T         
        %s: argument needed! All 9999 files like %sXXXX.fits exists! Brightness levels amount shoul be from 2 to 255 Can't add record to list Can't close mmap'ed file Can't copy data Can't get settings Can't make group operations 'in place'! Can't munmap Can't open %s for reading Can't read HDU Can't read file %s from sequence Can't read input file! Can't remove file %s Can't set settings Can't setup console Can't stat %s Conversion %s <%s> parameters:
 Degrees should be less than 360 Error in group operation Error on pipeline processing! FFTW error Filter FWHM should be not less than 1. Filter window height changed to 5 Filter window width changed to 5 Found %d pixels with undefined value Group operations need at least two FITS-files Image statistics after pipeline:
 Images with > 2 dimensions are not supported Input image statistics:
 Integer out of range Missed input file name! Mmap error for input No filename given! No input image given No pipeline parameters given Parameter '-i' isn't used in group operations! Pipeline parameters:
 Posterisation filter parameters should be:
%s Prewitt horizontal - simple derivative Prewitt vertical Scharr (modified Sobel) horizontal Scharr vertical Set output file name (-o) or its prefix (without key) The amount of available files less than two The output file exists Unknown group operation Usage: %s [args] [outfile prefix] [file list for group operations]

	Where args are:
 Wrong argument "%s" of parameter "%s" Wrong double number format! Wrong helpstring! Wrong levels amount: %d Wrong parameter: %s Wrong pipeline 'type' parameter: %s Wrong pipeline parameters! You should set 'scale' parameter add record to FITS-header arguments are absent calculate mean of specified FITS-files calculate median of specified FITS-files calculate sum of specified FITS-files delete all records with given substring delete given key gaussian horizontal Sobel input file laplasian of gaussian median filter nsteps	amount of steps
scale	scale type (uniform, log, exp, sqrt, pow) output file posterisation r	radius of filter (uint, 0 for cross 3x3) rewrite output file if exists (works only with option -i) save results into same file set pipeline parameters, arg: type=type:[help]:...
		type - transformation type (help for list)
		help - list of available parameters for given 'type' show statistic parameters of input and output image show this help sigma_x is too low, set to 1. sigma_y is too low, set to 1. simple adaptive median filter simple gradient (by Sobel) sx,sy	sigma by axes x & y
w,h	non-zero window width & height verbose level (each -v increase it) vertical Sobel Project-Id-Version: PACKAGE VERSION
Report-Msgid-Bugs-To: 
POT-Creation-Date: 2018-11-26 23:08+0300
PO-Revision-Date: YEAR-MO-DA HO:MI+ZONE
Last-Translator: FULL NAME <EMAIL@ADDRESS>
Language-Team: LANGUAGE <LL@li.org>
Language: 
MIME-Version: 1.0
Content-Type: text/plain; charset=koi8-r
Content-Transfer-Encoding: 8bit
 %s: ��������� ��������! ��� 9999 ������ ���� %sXXXX.fits ������! ���������� ������� �������� ������� ������ ���� �� 2 �� 255 �� ���� �������� ������ � ������ �� ���� ������� ��������� ���� �� ���� ����������� ������ ���������� �������� ��������� �� ���� ��������� ��������� �������� � ����������� '� ��� �� ����'! �� ���� ������� munmap �� ���� ��������� %s �� ���� �������� HDU �� ���� �������� ���� %s �� ������������������ ���������� �������� ������� ����! �� ���� ������� ��������� ���� �� ���� ��������� ��������� �� ���� ��������� �������, �� ���� �������� ��������� %s �������������� %s: %s
 ���� ������ ���� ������ 360 �������� ������ ��� ��������� ��������� ������������ ��������: %s ������ FFTW ���������� ������� ������ ���� �� ������ 1. ������ ������� �������� �� 5 ������ ������� �������� �� 5 ���������� %d �������� � �������������� ��������� ��������� �������� ������� ������� ��� ������� ���� FITS-������ ���������� �� ����������� ����� ���������:
 ����������� ����������� �� �������������� ���������� �� �������� �����������:
 ����� ��� ����������� ��������� �� ������ ��� �������� ����� ������ mmap �� ������ ��� ����� �� ������ ��� ����� ������������ ��������: %s �������� '-i' �� ������������ ��� ��������� ���������! ��������� ���������:
 ��������� ������� ������������ ������ ���� ������:
%s �������������� ������ ������� (���������� �����������) ������������ ������ ������� �������������� ������ ����� (���������������� ������) ������������ ������ ����� ������� ��� ��������� ����� (-o) ��� ��� ������� (��� �����) ��������� ��������� ������ ������ ���� �������� ���� ���������� ����������� ��������� �������� �������������: %s [���������] [������� ��������� �����]

	��� ���������:
 ������������ �������� "%s" ��������� "%s" ������������ ������ ����� double! ������������ ������ ������ ������ ������������ ���������� �������: %d ������������ ��������: %s ������������ ��������: %s ����������� ������ ��������� ��������� �� ������ �������� scale �������� ������ � FITS-����� ��������� ����������� ��������� ������� �������������� ���� ������������� ����������� ��������� ������� ���� ������������� ����������� ��������� ����� ���� ������������� ����������� ������� ��� ������ � ��������� ���������� ������� ��������� ���� ������� ������ �������������� ������ ������ ������� ���� ��������� ��������� ��������� ������ nsteps	���������� �������� �������
scale	������� �������������� (uniform, log, exp, sqrt, pow) �������� ���� "������������" r	������ ������� (����������� ����� ��� 0 ��� "������" 3x3) ������������ �������� ����, ���� �� ���������� (������ � ������ -i) �������� ��������� � ��� �� ���� ���������� ��������� ���������, ���������: type:[help]:...
		type - ��� �������������� (help ��� �������)
		help - ������ ��������� ��� ������� 'type' ����� ���������� �������������� ��������� �������� � ��������� ����������� ���������� ��� ������� sigma_x ������� ����, ������������ � 1. sigma_x ������� ����, ������������ � 1. ������� ���������� ��������� ������ ������� �������� (����������� ������) sx,sy	�������� ����� �� ���� x � y
w,h	������ � ������ ���������� ���� ������� ������� ������������ ������ (������ -v ����������� ���) ������������ ������ ������ 