function [ pdfText ] = pdfRead( pdf_location )
% PDFREAD reads a PDF file using the iText java library.
%
%   INPUT:
%       PDF_LOCATION:
%       String specifying the location of the PDF.
%
%   OUTPUT:
%       PDFTEXT:
%       Cell array, each cell corresponds to each page of the parsed PDF
%       file. Images are not extracted, only text.
%
% D. Wood, 07/03/2017

% add the java path for the iText java library (cedric)
jvpath = fileparts(which('iText-4.2.0-com.itextpdf.jar'));
javaaddpath(fullfile(jvpath,'iText-4.2.0-com.itextpdf.jar'));

% Start PDF reader object
reader = com.itextpdf.text.pdf.PdfReader(pdf_location);
tx = com.itextpdf.text.pdf.parser.PdfTextExtractor(reader);

% read each page, return output.
pdfText = cell(1);
page_num = 1;
page_limit_exceeded = false;
while ~page_limit_exceeded
    try
        t1 = tx.getTextFromPage(page_num);
        pdfText{page_num} = char(t1);
        page_num = page_num + 1;
    catch
        page_limit_exceeded = true;
    end
end
reader.close();



end

