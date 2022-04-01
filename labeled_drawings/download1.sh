#!/bin/bash

wget -P part1/ https://images.template.net/wp-content/uploads/2015/04/Cartoon-Lion.jpg
wget -P part1/ https://images.template.net/wp-content/uploads/2015/04/Cartoon-Pencil-Drawings.jpg
wget -P part1/ https://www.drawingforall.net/wp-content/uploads/2015/11/7-drawing-cartoons.jpg
wget -P part1/ https://www.easy-drawings-and-sketches.com/images/how-to-draw-cartoon-dog14.jpg
wget -P part1/ https://www.easy-drawings-and-sketches.com/images/cartoon-lion07.jpg
wget -P part1/ https://www.easy-drawings-and-sketches.com/images/cartoon-tiger06.jpg
wget -P part1/ https://www.easy-drawings-and-sketches.com/images/cartoon-bird06.jpg
wget -P part1/ https://www.easy-drawings-and-sketches.com/images/how-to-draw-cyclamen01.jpg
wget -P part1/ https://www.easy-drawings-and-sketches.com/images/how-to-draw-a-fish10.jpg
wget -P part1/ https://www.easy-drawings-and-sketches.com/images/how-to-draw-whale07.jpg
wget -P part1/ https://www.easy-drawings-and-sketches.com/images/how-to-draw-a-frog05.jpg
wget -P part1/ https://www.easy-drawings-and-sketches.com/images/draw-a-pig07.jpg
wget -P part1/ https://www.easy-drawings-and-sketches.com/images/draw-mickey-mouse07.jpg
wget -P part1/ https://www.easy-drawings-and-sketches.com/images/draw-minnie08.jpg
wget -P part1/ https://www.easy-drawings-and-sketches.com/images/how-to-draw-bambi08.jpg
wget -P part1/ https://www.easy-drawings-and-sketches.com/images/how-to-draw-dumbo07.jpg
wget -P part1/ https://www.easy-drawings-and-sketches.com/images/how-to-draw-nemo05.jpg
wget -P part1/ https://www.easy-drawings-and-sketches.com/images/how-to-draw-daisy-duck09.jpg
wget -P part1/ https://www.easy-drawings-and-sketches.com/images/how-to-draw-goofy07.jpg
wget -P part1/ https://www.easy-drawings-and-sketches.com/images/how-to-draw-donald-duck07.jpg
wget -P part1/ https://www.easy-drawings-and-sketches.com/images/cartoon-donkey07.jpg
wget -P part1/ https://www.allaboutdrawings.com/image-files/xcartoon-fish-drawing.jpg.pagespeed.ic.4HiG9-4vZJ.jpg
wget -P part1/ https://www.allaboutdrawings.com/image-files/cartoon-whale-drawing.jpg
wget -P part1/ https://www.easy-drawings-and-sketches.com/images/cartoon-pig09.jpg
wget -P part1/ https://www.easy-drawings-and-sketches.com/images/how-to-draw-popeye05.jpg
wget -P part1/ https://www.easy-drawings-and-sketches.com/images/how-to-draw-pikachu10.jpg
wget -P part1/ https://www.easy-drawings-and-sketches.com/images/cartoon-fish-drawing07.jpg
wget -P part1/ https://www.easy-drawings-and-sketches.com/images/how-to-draw-tinkerbell06.jpg
wget -P part1/ https://www.easy-drawings-and-sketches.com/images/how-to-draw-spongebob06.jpg
wget -P part1/ https://www.drawingforall.net/wp-content/uploads/2020/11/11-how-to-draw-a-cat-easy.jpg
wget -P part1/ https://www.ejin.ru/wp-content/uploads/2018/12/3d125deb8c18db9662098cbaeb6e7aa2-442x800.jpg
wget -P part1/ https://www.ejin.ru/wp-content/uploads/2018/10/5-5-618x800.jpg
wget -P part1/ https://www.ejin.ru/wp-content/uploads/2018/10/be_8-1-800x579.jpg
wget -P part1/ https://www.ejin.ru/wp-content/uploads/2018/10/orig-2-1.jpeg
wget -P part1/ https://www.ejin.ru/wp-content/uploads/2018/10/kak_narisovat_vinil_skretch_iz_my_little_pony_karandashom_pojetapno-7-600x800.jpg
wget -P part1/ https://www.ejin.ru/wp-content/uploads/2018/10/uchimsya_risovat_milogo_anime_lisenoka_prostym_karandashom-5-696x800.jpg
wget -P part1/ https://www.ejin.ru/wp-content/uploads/2018/10/5-33-800x588.jpg
wget -P part1/ https://www.ejin.ru/wp-content/uploads/2018/10/scrn_big_1-800x523.jpg
wget -P part1/ https://www.ejin.ru/wp-content/uploads/2018/10/ejsvq6z35b8-768x1094-562x800.jpg
wget -P part1/ https://www.ejin.ru/wp-content/uploads/2018/10/orig-1-5.jpg
wget -P part1/ https://www.ejin.ru/wp-content/uploads/2018/10/risunki-karandashom-na-novyj-god-2017-shablon-14-800x579.jpg
wget -P part1/ https://www.ejin.ru/wp-content/uploads/2018/08/kot-e1534886674525-696x435.jpg
wget -P part1/ https://www.ejin.ru/wp-content/uploads/2018/08/utenok_3797.jpg
wget -P part1/ https://www.ejin.ru/wp-content/uploads/2018/10/leonia_kitten-734x800.jpg
wget -P part1/ https://www.ejin.ru/wp-content/uploads/2018/10/d2eae9f0fc950993d81bd23d00365dad-e1506674135134-763x800.png
wget -P part1/ https://www.ejin.ru/wp-content/uploads/2018/10/neko_dance_by_eternalnova.jpg
wget -P part1/ https://www.ejin.ru/wp-content/uploads/2018/08/104246.jpg
wget -P part1/ https://www.ejin.ru/wp-content/uploads/2018/08/raskraski-graviti-folz17-800x647.jpg
wget -P part1/ https://www.ejin.ru/wp-content/uploads/2018/08/raskraski-iz-multikov-graviti-folz-14-600x800.jpg

mogrify -format png part1/*.jpg
mogrify -format png part1/*.jpeg

rm part1/*.jpg
rm part1/*.jpeg

a=1
for fff in part1/*.png; do
  new=$(printf "part1/%04d.png" "$a") #04 pad to length of 4
  mv -i -- "$fff" "$new"
  let a=a+1
done