% PDFREAD Demo
% D. Wood, 7/3/2017
%

% add the java path for the iText java library.
javaaddpath('iText-4.2.0-com.itextpdf.jar')

% read the PDF file.
pdf_text = pdfRead('sample.pdf');

% the output is a cell array.
% each cell contains the text from a different page of the PDF file.
disp(pdf_text{1})


