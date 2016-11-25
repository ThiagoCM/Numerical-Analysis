%% SVD COMPRESSION %%
img = imread('_MG_7228.cr2'); %Lendo a imagem raw e salvando a mesma em img
figure(1); imshow(img, []);
title('Imagem Original');

%separando a imagem em três novas matrizes que representam as colorações
%rgb
red = double(img(:, :, 1));
green = double(img(:, :, 2));
blue = double(img(:, :, 3));

%plot das imagens separadas por sua coloração
figure(2);
subplot(2, 2, 1); imshow(red, []);
subplot(2, 2, 2); imshow(green, []);
subplot(2, 2, 3); imshow(blue, []);


%comprimindo pela técnica SVD cada uma das matrizes que representam
%vermelho, verde e azul da imagem original
[ur36, sr36, vr36] = svds(red, 36);
img_red36 = uint8(ur36 * sr36 * transpose(vr36));

[ug36, sg36, vg36] = svds(green, 36);
img_green36 = uint8(ug36 * sg36 * transpose(vg36));

[ub36, sb36, vb36] = svds(blue, 36);
img_blue36 = uint8(ub36 * sb36 * transpose(vb36));

img_com36(:, :, 1) = img_red36;
img_com36(:, :, 2) = img_green36;
img_com36(:, :, 3) = img_blue36;

figure(3);
subplot(2, 2, 1); imshow(img_red36, []);
subplot(2, 2, 2); imshow(img_green36, []);
subplot(2, 2, 3); imshow(img_blue36, []);

figure(4);
imshow(img_com36, []);imsave;
title('Imagem Comprimida, 36 maiores numeros');

[ur48, sr48, vr48] = svds(red, 48);
img_red48 = uint8(ur48 * sr48 * transpose(vr48));

[ug48, sg48, vg48] = svds(green, 48);
img_green48 = uint8(ug48 * sg48 * transpose(vg48));

[ub48, sb48, vb48] = svds(blue, 48);
img_blue48 = uint8(ub48 * sb48 * transpose(vb48));

img_com48(:, :, 1) = img_red48;
img_com48(:, :, 2) = img_green48;
img_com48(:, :, 3) = img_blue48;

figure(5);
subplot(2, 2, 1); imshow(img_red48, []);
subplot(2, 2, 2); imshow(img_green48, []);
subplot(2, 2, 3); imshow(img_blue48, []);

figure(6);
imshow(img_com48, []);imsave;
title('Imagem Comprimida, 48 maiores numeros');
