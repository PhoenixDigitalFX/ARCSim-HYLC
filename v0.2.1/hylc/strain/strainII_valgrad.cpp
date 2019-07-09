#include "strain.hpp"
#ifdef hylc_strain_II

using namespace hylc;
using namespace hylc::mathematica;

std::tuple<Mat6x18, Vec6> hylc::mathematica::strain_valgrad(
    const Vec18 &xloc, const Mat2x2 &invDm, const Real &A,
    const Real &thetarest0, const Real &thetarest1, const Real &thetarest2,
    const Real &l0, const Real &l1, const Real &l2, const Vec2 &t0,
    const Vec2 &t1, const Vec2 &t2) {
  // define output
  Mat6x18 grad(0);
  Vec6 val(0);
  auto out1 = [&](int i) -> Real & { return val[i]; };
  auto out2 = [&](int i, int j) -> Real & { return grad(i, j); };

  Real copt57  = invDm(0, 0);
  Real copt162 = xloc(0);
  Real copt232 = -copt162;
  Real copt300 = xloc(3);
  Real copt301 = copt232 + copt300;
  Real copt302 = copt301 * copt57;
  Real copt304 = invDm(1, 0);
  Real copt306 = xloc(6);
  Real copt308 = copt232 + copt306;
  Real copt309 = copt304 * copt308;
  Real copt310 = copt302 + copt309;
  Real copt311 = Power(copt310, 2);
  Real copt317 = xloc(1);
  Real copt399 = -copt317;
  Real copt410 = xloc(4);
  Real copt416 = copt399 + copt410;
  Real copt417 = copt416 * copt57;
  Real copt418 = xloc(7);
  Real copt438 = copt399 + copt418;
  Real copt439 = copt304 * copt438;
  Real copt440 = copt417 + copt439;
  Real copt442 = Power(copt440, 2);
  Real copt443 = xloc(2);
  Real copt445 = -copt443;
  Real copt446 = xloc(5);
  Real copt447 = copt445 + copt446;
  Real copt448 = copt447 * copt57;
  Real copt449 = xloc(8);
  Real copt450 = copt445 + copt449;
  Real copt451 = copt304 * copt450;
  Real copt452 = copt448 + copt451;
  Real copt453 = Power(copt452, 2);
  Real copt454 = copt311 + copt442 + copt453;
  Real copt455 = Sqrt(copt454);
  Real copt456 = 1 / copt455;
  Real copt457 = invDm(0, 1);
  Real copt459 = invDm(1, 1);
  Real copt458 = copt301 * copt457;
  Real copt460 = copt308 * copt459;
  Real copt461 = copt458 + copt460;
  Real copt472 = Power(copt461, 2);
  Real copt463 = copt416 * copt457;
  Real copt464 = copt438 * copt459;
  Real copt465 = copt463 + copt464;
  Real copt473 = Power(copt465, 2);
  Real copt467 = copt447 * copt457;
  Real copt468 = copt450 * copt459;
  Real copt469 = copt467 + copt468;
  Real copt474 = Power(copt469, 2);
  Real copt475 = copt472 + copt473 + copt474;
  Real copt476 = Sqrt(copt475);
  Real copt477 = 1 / copt476;
  Real copt479 = 1 / A;
  Real copt480 = 1 / l0;
  Real copt481 = 1 / l1;
  Real copt482 = 1 / l2;
  Real copt492 = xloc(9);
  Real copt497 = xloc(10);
  Real copt484 = -copt300;
  Real copt516 = xloc(11);
  Real copt499 = -copt410;
  Real copt509 = -copt446;
  Real copt541 = copt162 + copt484;
  Real copt542 = Power(copt541, 2);
  Real copt543 = copt317 + copt499;
  Real copt544 = Power(copt543, 2);
  Real copt545 = copt443 + copt509;
  Real copt546 = Power(copt545, 2);
  Real copt547 = copt542 + copt544 + copt546;
  Real copt548 = Sqrt(copt547);
  Real copt504 = -copt306;
  Real copt531 = -copt497;
  Real copt493 = -copt492;
  Real copt526 = -copt449;
  Real copt505 = copt300 + copt504;
  Real copt524 = copt418 + copt499;
  Real copt589 = xloc(12);
  Real copt593 = xloc(13);
  Real copt506 = copt443 * copt505;
  Real copt507 = copt306 * copt446;
  Real copt508 = -(copt300 * copt449);
  Real copt510 = copt449 + copt509;
  Real copt511 = copt162 * copt510;
  Real copt512 = copt506 + copt507 + copt508 + copt511;
  Real copt591 = copt504 + copt589;
  Real copt602 = xloc(14);
  Real copt487 = -copt418;
  Real copt488 = copt410 + copt487;
  Real copt603 = -copt602;
  Real copt604 = copt449 + copt603;
  Real copt622 = Power(copt505, 2);
  Real copt623 = Power(copt488, 2);
  Real copt527 = copt446 + copt526;
  Real copt624 = Power(copt527, 2);
  Real copt625 = copt622 + copt623 + copt624;
  Real copt626 = Sqrt(copt625);
  Real copt549 = -(copt162 * copt418 * copt446);
  Real copt550 = copt162 * copt410 * copt449;
  Real copt590 = -(copt418 * copt589);
  Real copt592 = copt410 * copt591;
  Real copt594 = -copt593;
  Real copt595 = copt418 + copt594;
  Real copt596 = copt300 * copt595;
  Real copt597 = copt306 * copt593;
  Real copt598 = copt590 + copt592 + copt596 + copt597;
  Real copt584 = copt317 * copt505;
  Real copt585 = copt306 * copt410;
  Real copt586 = -(copt300 * copt418);
  Real copt587 = copt162 * copt524;
  Real copt588 = copt584 + copt585 + copt586 + copt587;
  Real copt654 = xloc(15);
  Real copt659 = xloc(16);
  Real copt655 = -copt654;
  Real copt656 = copt306 + copt655;
  Real copt667 = xloc(17);
  Real copt609 = copt443 * copt488;
  Real copt610 = copt418 * copt446;
  Real copt611 = -(copt410 * copt449);
  Real copt612 = copt317 * copt510;
  Real copt613 = copt609 + copt610 + copt611 + copt612;
  Real copt669 = copt526 + copt667;
  Real copt682 = copt162 + copt504;
  Real copt683 = Power(copt682, 2);
  Real copt684 = copt317 + copt487;
  Real copt685 = Power(copt684, 2);
  Real copt686 = copt443 + copt526;
  Real copt687 = Power(copt686, 2);
  Real copt688 = copt683 + copt685 + copt687;
  Real copt689 = Sqrt(copt688);
  Real copt658 = copt418 * copt654;
  Real copt660 = -(copt306 * copt659);
  Real copt661 = copt487 + copt659;
  Real copt462 = copt310 * copt461;
  Real copt466 = copt440 * copt465;
  Real copt470 = copt452 * copt469;
  Real copt471 = copt462 + copt466 + copt470;
  Real copt483 = -(copt306 * copt410);
  Real copt485 = copt306 + copt484;
  Real copt486 = copt317 * copt485;
  Real copt489 = copt162 * copt488;
  Real copt490 = copt300 * copt418;
  Real copt491 = copt483 + copt486 + copt489 + copt490;
  Real copt494 = copt300 + copt493;
  Real copt495 = copt317 * copt494;
  Real copt496 = copt410 * copt492;
  Real copt498 = -(copt300 * copt497);
  Real copt500 = copt497 + copt499;
  Real copt501 = copt162 * copt500;
  Real copt502 = copt495 + copt496 + copt498 + copt501;
  Real copt503 = copt491 * copt502;
  Real copt513 = -(copt446 * copt492);
  Real copt514 = copt484 + copt492;
  Real copt515 = copt443 * copt514;
  Real copt517 = -copt516;
  Real copt518 = copt446 + copt517;
  Real copt519 = copt162 * copt518;
  Real copt520 = copt300 * copt516;
  Real copt521 = copt513 + copt515 + copt519 + copt520;
  Real copt522 = copt512 * copt521;
  Real copt523 = -(copt418 * copt446);
  Real copt525 = copt443 * copt524;
  Real copt528 = copt317 * copt527;
  Real copt529 = copt410 * copt449;
  Real copt530 = copt523 + copt525 + copt528 + copt529;
  Real copt532 = copt410 + copt531;
  Real copt533 = copt443 * copt532;
  Real copt534 = copt446 * copt497;
  Real copt535 = -(copt410 * copt516);
  Real copt536 = copt509 + copt516;
  Real copt537 = copt317 * copt536;
  Real copt538 = copt533 + copt534 + copt535 + copt537;
  Real copt539 = copt530 * copt538;
  Real copt540 = copt503 + copt522 + copt539;
  Real copt551 = copt418 * copt446 * copt492;
  Real copt552 = -(copt410 * copt449 * copt492);
  Real copt553 = copt162 * copt446 * copt497;
  Real copt554 = -(copt306 * copt446 * copt497);
  Real copt555 = -(copt162 * copt449 * copt497);
  Real copt556 = copt300 * copt449 * copt497;
  Real copt557 = -(copt418 * copt492);
  Real copt558 = copt492 + copt504;
  Real copt559 = copt410 * copt558;
  Real copt560 = copt418 + copt531;
  Real copt561 = copt300 * copt560;
  Real copt562 = copt306 * copt497;
  Real copt563 = copt557 + copt559 + copt561 + copt562;
  Real copt564 = copt443 * copt563;
  Real copt565 = -(copt162 * copt410 * copt516);
  Real copt566 = copt306 * copt410 * copt516;
  Real copt567 = copt162 * copt418 * copt516;
  Real copt568 = -(copt300 * copt418 * copt516);
  Real copt569 = copt306 + copt493;
  Real copt570 = copt446 * copt569;
  Real copt571 = copt449 * copt492;
  Real copt572 = -(copt306 * copt516);
  Real copt573 = copt516 + copt526;
  Real copt574 = copt300 * copt573;
  Real copt575 = copt570 + copt571 + copt572 + copt574;
  Real copt576 = copt317 * copt575;
  Real copt577 = copt549 + copt550 + copt551 + copt552 + copt553 + copt554 +
                 copt555 + copt556 + copt564 + copt565 + copt566 + copt567 +
                 copt568 + copt576;
  Real copt578 = copt548 * copt577;
  Real copt579 = ArcTan(copt540, copt578);
  Real copt580 = copt579 + thetarest0;
  Real copt581 = t0(0);
  Real copt724 = Power(copt581, 2);
  Real copt599 = copt588 * copt598;
  Real copt600 = -(copt449 * copt589);
  Real copt601 = copt446 * copt591;
  Real copt605 = copt300 * copt604;
  Real copt606 = copt306 * copt602;
  Real copt607 = copt600 + copt601 + copt605 + copt606;
  Real copt608 = copt512 * copt607;
  Real copt614 = -(copt449 * copt593);
  Real copt615 = copt487 + copt593;
  Real copt616 = copt446 * copt615;
  Real copt617 = copt410 * copt604;
  Real copt618 = copt418 * copt602;
  Real copt619 = copt614 + copt616 + copt617 + copt618;
  Real copt620 = copt613 * copt619;
  Real copt621 = copt599 + copt608 + copt620;
  Real copt627 = copt418 * copt446 * copt589;
  Real copt628 = -(copt410 * copt449 * copt589);
  Real copt629 = copt162 * copt446 * copt593;
  Real copt630 = -(copt306 * copt446 * copt593);
  Real copt631 = -(copt162 * copt449 * copt593);
  Real copt632 = copt300 * copt449 * copt593;
  Real copt633 = copt443 * copt598;
  Real copt634 = -(copt162 * copt410 * copt602);
  Real copt635 = copt306 * copt410 * copt602;
  Real copt636 = copt162 * copt418 * copt602;
  Real copt637 = -(copt300 * copt418 * copt602);
  Real copt638 = -copt589;
  Real copt639 = copt306 + copt638;
  Real copt640 = copt446 * copt639;
  Real copt641 = copt449 * copt589;
  Real copt642 = -(copt306 * copt602);
  Real copt643 = copt526 + copt602;
  Real copt644 = copt300 * copt643;
  Real copt645 = copt640 + copt641 + copt642 + copt644;
  Real copt646 = copt317 * copt645;
  Real copt647 = copt549 + copt550 + copt627 + copt628 + copt629 + copt630 +
                 copt631 + copt632 + copt633 + copt634 + copt635 + copt636 +
                 copt637 + copt646;
  Real copt648 = copt626 * copt647;
  Real copt649 = ArcTan(copt621, copt648);
  Real copt650 = copt649 + thetarest1;
  Real copt651 = t1(0);
  Real copt726 = Power(copt651, 2);
  Real copt657 = copt317 * copt656;
  Real copt662 = copt162 * copt661;
  Real copt663 = copt657 + copt658 + copt660 + copt662;
  Real copt664 = copt588 * copt663;
  Real copt665 = copt443 * copt656;
  Real copt666 = copt449 * copt654;
  Real copt668 = -(copt306 * copt667);
  Real copt670 = copt162 * copt669;
  Real copt671 = copt665 + copt666 + copt668 + copt670;
  Real copt672 = copt512 * copt671;
  Real copt673 = -copt659;
  Real copt674 = copt418 + copt673;
  Real copt675 = copt443 * copt674;
  Real copt676 = copt449 * copt659;
  Real copt677 = -(copt418 * copt667);
  Real copt678 = copt317 * copt669;
  Real copt679 = copt675 + copt676 + copt677 + copt678;
  Real copt680 = copt613 * copt679;
  Real copt681 = copt664 + copt672 + copt680;
  Real copt690 = copt162 * copt418 * copt446;
  Real copt691 = -(copt162 * copt410 * copt449);
  Real copt692 = -(copt418 * copt446 * copt654);
  Real copt693 = copt410 * copt449 * copt654;
  Real copt694 = -(copt162 * copt446 * copt659);
  Real copt695 = copt306 * copt446 * copt659;
  Real copt696 = copt162 * copt449 * copt659;
  Real copt697 = -(copt300 * copt449 * copt659);
  Real copt698 = copt410 * copt656;
  Real copt699 = copt300 * copt661;
  Real copt700 = copt658 + copt660 + copt698 + copt699;
  Real copt701 = copt443 * copt700;
  Real copt702 = copt162 * copt410 * copt667;
  Real copt703 = -(copt306 * copt410 * copt667);
  Real copt704 = -(copt162 * copt418 * copt667);
  Real copt705 = copt300 * copt418 * copt667;
  Real copt706 = -(copt449 * copt654);
  Real copt707 = copt504 + copt654;
  Real copt708 = copt446 * copt707;
  Real copt709 = -copt667;
  Real copt710 = copt449 + copt709;
  Real copt711 = copt300 * copt710;
  Real copt712 = copt306 * copt667;
  Real copt713 = copt706 + copt708 + copt711 + copt712;
  Real copt714 = copt317 * copt713;
  Real copt715 = copt690 + copt691 + copt692 + copt693 + copt694 + copt695 +
                 copt696 + copt697 + copt701 + copt702 + copt703 + copt704 +
                 copt705 + copt714;
  Real copt716 = -(copt689 * copt715);
  Real copt717 = ArcTan(copt681, copt716);
  Real copt718 = copt717 + thetarest2;
  Real copt719 = t2(0);
  Real copt728 = Power(copt719, 2);
  Real copt733 = Power(copt471, 2);
  Real copt734 = -copt733;
  Real copt735 = copt454 * copt475;
  Real copt736 = copt734 + copt735;
  Real copt737 = 1 / copt736;
  Real copt739 = copt304 * copt457;
  Real copt740 = -(copt459 * copt57);
  Real copt741 = copt739 + copt740;
  Real copt742 = Power(copt741, 2);
  Real copt743 = 1 / copt742;
  Real copt744 = Power(copt162, 2);
  Real copt745 = Power(copt410, 2);
  Real copt746 = copt744 * copt745;
  Real copt747 = Power(copt446, 2);
  Real copt748 = copt744 * copt747;
  Real copt749 = -2 * copt162 * copt306 * copt745;
  Real copt750 = -2 * copt162 * copt306 * copt747;
  Real copt751 = Power(copt306, 2);
  Real copt752 = copt745 * copt751;
  Real copt753 = copt747 * copt751;
  Real copt754 = -2 * copt410 * copt418 * copt744;
  Real copt755 = 2 * copt162 * copt300 * copt410 * copt418;
  Real copt756 = 2 * copt162 * copt306 * copt410 * copt418;
  Real copt757 = -2 * copt300 * copt306 * copt410 * copt418;
  Real copt758 = Power(copt418, 2);
  Real copt759 = copt744 * copt758;
  Real copt760 = -2 * copt162 * copt300 * copt758;
  Real copt761 = Power(copt300, 2);
  Real copt762 = copt758 * copt761;
  Real copt763 = copt747 * copt758;
  Real copt764 = Power(copt443, 2);
  Real copt765 = -2 * copt300 * copt306;
  Real copt766 = -2 * copt410 * copt418;
  Real copt767 = copt745 + copt751 + copt758 + copt761 + copt765 + copt766;
  Real copt768 = copt764 * copt767;
  Real copt769 = -2 * copt446 * copt449 * copt744;
  Real copt770 = 2 * copt162 * copt300 * copt446 * copt449;
  Real copt771 = 2 * copt162 * copt306 * copt446 * copt449;
  Real copt772 = -2 * copt300 * copt306 * copt446 * copt449;
  Real copt773 = -2 * copt410 * copt418 * copt446 * copt449;
  Real copt774 = Power(copt449, 2);
  Real copt775 = copt744 * copt774;
  Real copt776 = -2 * copt162 * copt300 * copt774;
  Real copt777 = copt761 * copt774;
  Real copt778 = copt745 * copt774;
  Real copt779 = Power(copt317, 2);
  Real copt780 = -2 * copt446 * copt449;
  Real copt781 = copt747 + copt751 + copt761 + copt765 + copt774 + copt780;
  Real copt782 = copt779 * copt781;
  Real copt783 = -(copt300 * copt306 * copt410);
  Real copt784 = copt410 * copt751;
  Real copt785 = copt162 * copt488 * copt505;
  Real copt786 = copt418 * copt761;
  Real copt787 = copt418 * copt747;
  Real copt788 = -(copt300 * copt306 * copt418);
  Real copt789 = copt443 * copt488 * copt527;
  Real copt790 = -(copt410 * copt446 * copt449);
  Real copt791 = -(copt418 * copt446 * copt449);
  Real copt792 = copt410 * copt774;
  Real copt793 = copt783 + copt784 + copt785 + copt786 + copt787 + copt788 +
                 copt789 + copt790 + copt791 + copt792;
  Real copt794 = -2 * copt317 * copt793;
  Real copt795 = copt446 * copt751;
  Real copt796 = -(copt410 * copt418 * copt446);
  Real copt797 = copt446 * copt758;
  Real copt798 = copt162 * copt505 * copt527;
  Real copt799 = copt449 * copt761;
  Real copt800 = copt449 * copt745;
  Real copt801 = -(copt410 * copt418 * copt449);
  Real copt802 = copt446 + copt449;
  Real copt803 = -(copt300 * copt306 * copt802);
  Real copt804 = copt795 + copt796 + copt797 + copt798 + copt799 + copt800 +
                 copt801 + copt803;
  Real copt805 = -2 * copt443 * copt804;
  Real copt806 = copt746 + copt748 + copt749 + copt750 + copt752 + copt753 +
                 copt754 + copt755 + copt756 + copt757 + copt759 + copt760 +
                 copt762 + copt763 + copt768 + copt769 + copt770 + copt771 +
                 copt772 + copt773 + copt775 + copt776 + copt777 + copt778 +
                 copt782 + copt794 + copt805;
  Real copt807 = 1 / copt806;
  Real copt582 = t0(1);
  Real copt652 = t1(1);
  Real copt720 = t2(1);
  Real copt815 = Power(copt457, 2);
  Real copt832 = Power(copt459, 2);
  Real copt840 = Power(copt582, 2);
  Real copt843 = Power(copt652, 2);
  Real copt846 = Power(copt720, 2);
  Real copt816 = -2 * copt162 * copt300;
  Real copt817 = -2 * copt317 * copt410;
  Real copt818 = -2 * copt443 * copt446;
  Real copt819 = copt744 + copt745 + copt747 + copt761 + copt764 + copt779 +
                 copt816 + copt817 + copt818;
  Real copt821 = -(copt443 * copt446);
  Real copt822 = copt300 * copt306;
  Real copt823 = copt300 + copt306;
  Real copt824 = -(copt162 * copt823);
  Real copt825 = copt410 * copt418;
  Real copt826 = copt410 + copt418;
  Real copt827 = -(copt317 * copt826);
  Real copt828 = -(copt443 * copt449);
  Real copt829 = copt446 * copt449;
  Real copt830 = copt744 + copt764 + copt779 + copt821 + copt822 + copt824 +
                 copt825 + copt827 + copt828 + copt829;
  Real copt833 = -2 * copt162 * copt306;
  Real copt834 = -2 * copt317 * copt418;
  Real copt835 = -2 * copt443 * copt449;
  Real copt836 = copt744 + copt751 + copt758 + copt764 + copt774 + copt779 +
                 copt833 + copt834 + copt835;
  Real copt841  = copt840 * l1 * l2 * thetarest0;
  Real copt842  = copt579 * copt840 * l1 * l2;
  Real copt844  = copt843 * l0 * l2 * thetarest1;
  Real copt845  = copt649 * copt843 * l0 * l2;
  Real copt847  = copt846 * l0 * l1 * thetarest2;
  Real copt848  = copt717 * copt846 * l0 * l1;
  Real copt849  = copt841 + copt842 + copt844 + copt845 + copt847 + copt848;
  Real copt862  = Power(copt57, 2);
  Real copt865  = Power(copt304, 2);
  Real copt808  = copt581 * copt582 * l1 * l2 * thetarest0;
  Real copt809  = copt579 * copt581 * copt582 * l1 * l2;
  Real copt810  = copt651 * copt652 * l0 * l2 * thetarest1;
  Real copt811  = copt649 * copt651 * copt652 * l0 * l2;
  Real copt812  = copt719 * copt720 * l0 * l1 * thetarest2;
  Real copt813  = copt717 * copt719 * copt720 * l0 * l1;
  Real copt814  = copt808 + copt809 + copt810 + copt811 + copt812 + copt813;
  Real copt850  = copt457 * copt819;
  Real copt851  = copt459 * copt830;
  Real copt852  = copt850 + copt851;
  Real copt853  = copt57 * copt852;
  Real copt854  = copt457 * copt830;
  Real copt855  = copt459 * copt836;
  Real copt856  = copt854 + copt855;
  Real copt857  = copt304 * copt856;
  Real copt858  = copt853 + copt857;
  Real copt872  = copt304 + copt57;
  Real copt891  = copt454 * copt455;
  Real copt892  = 1 / copt891;
  Real copt893  = copt475 * copt476;
  Real copt894  = 1 / copt893;
  Real copt895  = copt457 + copt459;
  Real copt873  = copt541 * copt57;
  Real copt874  = copt304 * copt682;
  Real copt875  = copt873 + copt874;
  Real copt877  = copt543 * copt57;
  Real copt878  = copt304 * copt684;
  Real copt879  = copt877 + copt878;
  Real copt881  = copt545 * copt57;
  Real copt882  = copt304 * copt686;
  Real copt883  = copt881 + copt882;
  Real copt983  = -copt457;
  Real copt984  = -copt459;
  Real copt985  = copt983 + copt984;
  Real copt583  = copt580 * copt581 * copt582 * l1 * l2;
  Real copt653  = copt650 * copt651 * copt652 * l0 * l2;
  Real copt721  = copt718 * copt719 * copt720 * l0 * l1;
  Real copt722  = copt583 + copt653 + copt721;
  Real copt723  = copt471 * copt722;
  Real copt725  = copt580 * copt724 * l1 * l2;
  Real copt727  = copt650 * copt726 * l0 * l2;
  Real copt729  = copt718 * copt728 * l0 * l1;
  Real copt730  = copt725 + copt727 + copt729;
  Real copt731  = -(copt475 * copt730);
  Real copt732  = copt723 + copt731;
  Real copt897  = copt875 * copt895;
  Real copt898  = copt457 * copt541;
  Real copt899  = copt459 * copt682;
  Real copt900  = copt898 + copt899;
  Real copt901  = copt872 * copt900;
  Real copt902  = copt897 + copt901;
  Real copt1002 = Power(copt736, 2);
  Real copt1003 = 1 / copt1002;
  Real copt1012 = Power(copt540, 2);
  Real copt1013 = Power(copt577, 2);
  Real copt1014 = copt1013 * copt547;
  Real copt1015 = copt1012 + copt1014;
  Real copt1016 = 1 / copt1015;
  Real copt1022 = 1 / copt548;
  Real copt1031 = Power(copt621, 2);
  Real copt1032 = Power(copt647, 2);
  Real copt1033 = copt1032 * copt625;
  Real copt1034 = copt1031 + copt1033;
  Real copt1035 = 1 / copt1034;
  Real copt1048 = Power(copt715, 2);
  Real copt1049 = copt1048 * copt688;
  Real copt1050 = Power(copt681, 2);
  Real copt1051 = copt1049 + copt1050;
  Real copt1052 = 1 / copt1051;
  Real copt1058 = 1 / copt689;
  Real copt1007 = copt491 * copt500;
  Real copt1008 = copt488 * copt502;
  Real copt1009 = copt512 * copt518;
  Real copt1010 = copt510 * copt521;
  Real copt1011 = copt1007 + copt1008 + copt1009 + copt1010;
  Real copt1017 = -(copt1011 * copt1016 * copt548 * copt577);
  Real copt1018 = -(copt449 * copt497);
  Real copt1019 = copt418 * copt516;
  Real copt1020 = copt1018 + copt1019 + copt523 + copt529 + copt534 + copt535;
  Real copt1021 = copt1020 * copt548;
  Real copt1023 = copt1022 * copt541 * copt577;
  Real copt1024 = copt1021 + copt1023;
  Real copt1025 = copt1016 * copt1024 * copt540;
  Real copt1026 = copt1017 + copt1025;
  Real copt1028 = copt446 * copt593;
  Real copt1029 = -(copt410 * copt602);
  Real copt1030 = copt1028 + copt1029 + copt523 + copt529 + copt614 + copt618;
  Real copt1036 = copt1030 * copt1035 * copt621 * copt626;
  Real copt1037 = copt524 * copt598;
  Real copt1038 = copt510 * copt607;
  Real copt1039 = copt1037 + copt1038;
  Real copt1040 = -(copt1035 * copt1039 * copt626 * copt647);
  Real copt1041 = copt1036 + copt1040;
  Real copt1043 = copt588 * copt661;
  Real copt1044 = copt524 * copt663;
  Real copt1045 = copt512 * copt669;
  Real copt1046 = copt510 * copt671;
  Real copt1047 = copt1043 + copt1044 + copt1045 + copt1046;
  Real copt1053 = copt1047 * copt1052 * copt689 * copt715;
  Real copt1054 = -(copt446 * copt659);
  Real copt1055 = copt410 * copt667;
  Real copt1056 = copt1054 + copt1055 + copt610 + copt611 + copt676 + copt677;
  Real copt1057 = -(copt1056 * copt689);
  Real copt1059 = -(copt1058 * copt682 * copt715);
  Real copt1060 = copt1057 + copt1059;
  Real copt1061 = copt1052 * copt1060 * copt681;
  Real copt1062 = copt1053 + copt1061;
  Real copt908  = copt879 * copt895;
  Real copt909  = copt457 * copt543;
  Real copt910  = copt459 * copt684;
  Real copt911  = copt909 + copt910;
  Real copt912  = copt872 * copt911;
  Real copt913  = copt908 + copt912;
  Real copt997  = -copt57;
  Real copt998  = -copt304;
  Real copt999  = copt997 + copt998;
  Real copt1081 = copt491 * copt494;
  Real copt1082 = copt485 * copt502;
  Real copt1083 = copt530 * copt536;
  Real copt1084 = copt527 * copt538;
  Real copt1085 = copt1081 + copt1082 + copt1083 + copt1084;
  Real copt1086 = -(copt1016 * copt1085 * copt548 * copt577);
  Real copt1087 = copt548 * copt575;
  Real copt1088 = copt1022 * copt543 * copt577;
  Real copt1089 = copt1087 + copt1088;
  Real copt1090 = copt1016 * copt1089 * copt540;
  Real copt1091 = copt1086 + copt1090;
  Real copt1093 = copt1035 * copt621 * copt626 * copt645;
  Real copt1094 = copt505 * copt598;
  Real copt1095 = copt510 * copt619;
  Real copt1096 = copt1094 + copt1095;
  Real copt1097 = -(copt1035 * copt1096 * copt626 * copt647);
  Real copt1098 = copt1093 + copt1097;
  Real copt1100 = copt588 * copt656;
  Real copt1101 = copt505 * copt663;
  Real copt1102 = copt613 * copt669;
  Real copt1103 = copt510 * copt679;
  Real copt1104 = copt1100 + copt1101 + copt1102 + copt1103;
  Real copt1105 = copt1052 * copt1104 * copt689 * copt715;
  Real copt1106 = -(copt689 * copt713);
  Real copt1107 = -(copt1058 * copt684 * copt715);
  Real copt1108 = copt1106 + copt1107;
  Real copt1109 = copt1052 * copt1108 * copt681;
  Real copt1110 = copt1105 + copt1109;
  Real copt919  = copt883 * copt895;
  Real copt920  = copt457 * copt545;
  Real copt921  = copt459 * copt686;
  Real copt922  = copt920 + copt921;
  Real copt923  = copt872 * copt922;
  Real copt924  = copt919 + copt923;
  Real copt1129 = copt512 * copt514;
  Real copt1130 = copt530 * copt532;
  Real copt1131 = copt505 * copt521;
  Real copt1132 = copt524 * copt538;
  Real copt1133 = copt1129 + copt1130 + copt1131 + copt1132;
  Real copt1134 = -(copt1016 * copt1133 * copt548 * copt577);
  Real copt1135 = copt548 * copt563;
  Real copt1136 = copt1022 * copt545 * copt577;
  Real copt1137 = copt1135 + copt1136;
  Real copt1138 = copt1016 * copt1137 * copt540;
  Real copt1139 = copt1134 + copt1138;
  Real copt1141 = copt1035 * copt598 * copt621 * copt626;
  Real copt1142 = copt505 * copt607;
  Real copt1143 = copt488 * copt619;
  Real copt1144 = copt1142 + copt1143;
  Real copt1145 = -(copt1035 * copt1144 * copt626 * copt647);
  Real copt1146 = copt1141 + copt1145;
  Real copt1148 = copt512 * copt656;
  Real copt1149 = copt613 * copt674;
  Real copt1150 = copt505 * copt671;
  Real copt1151 = copt488 * copt679;
  Real copt1152 = copt1148 + copt1149 + copt1150 + copt1151;
  Real copt1153 = copt1052 * copt1152 * copt689 * copt715;
  Real copt1154 = -(copt689 * copt700);
  Real copt1155 = -(copt1058 * copt686 * copt715);
  Real copt1156 = copt1154 + copt1155;
  Real copt1157 = copt1052 * copt1156 * copt681;
  Real copt1158 = copt1153 + copt1157;
  Real copt930  = copt304 * copt308 * copt457;
  Real copt1171 = 2 * copt301 * copt457 * copt57;
  Real copt1172 = copt308 * copt459 * copt57;
  Real copt1173 = copt1171 + copt1172 + copt930;
  Real copt1211 = 1 / copt626;
  Real copt1180 = copt317 + copt531;
  Real copt1181 = copt1180 * copt491;
  Real copt1182 = copt438 * copt502;
  Real copt1183 = copt445 + copt516;
  Real copt1184 = copt1183 * copt512;
  Real copt1185 = copt521 * copt686;
  Real copt1186 = copt1181 + copt1182 + copt1184 + copt1185;
  Real copt1187 = -(copt1016 * copt1186 * copt548 * copt577);
  Real copt1188 = copt443 * copt560;
  Real copt1189 = copt449 * copt497;
  Real copt1190 = -(copt418 * copt516);
  Real copt1191 = copt317 * copt573;
  Real copt1192 = copt1188 + copt1189 + copt1190 + copt1191;
  Real copt1193 = copt1192 * copt548;
  Real copt1194 = -(copt1022 * copt541 * copt577);
  Real copt1195 = copt1193 + copt1194;
  Real copt1196 = copt1016 * copt1195 * copt540;
  Real copt1197 = copt1187 + copt1196;
  Real copt1199 = copt588 * copt595;
  Real copt1200 = copt598 * copt684;
  Real copt1201 = copt512 * copt604;
  Real copt1202 = copt607 * copt686;
  Real copt1203 = copt1199 + copt1200 + copt1201 + copt1202;
  Real copt1204 = -(copt1035 * copt1203 * copt626 * copt647);
  Real copt1205 = copt443 * copt595;
  Real copt1206 = copt449 * copt593;
  Real copt1207 = -(copt418 * copt602);
  Real copt1208 = copt317 * copt643;
  Real copt1209 = copt1205 + copt1206 + copt1207 + copt1208;
  Real copt1210 = copt1209 * copt626;
  Real copt1212 = copt1211 * copt505 * copt647;
  Real copt1213 = copt1210 + copt1212;
  Real copt1214 = copt1035 * copt1213 * copt621;
  Real copt1215 = copt1204 + copt1214;
  Real copt1217 = copt663 * copt684;
  Real copt1218 = copt671 * copt686;
  Real copt1219 = copt1217 + copt1218;
  Real copt1220 = copt1052 * copt1219 * copt689 * copt715;
  Real copt1221 = -(copt449 * copt659);
  Real copt1222 = copt443 * copt661;
  Real copt1223 = copt317 * copt710;
  Real copt1224 = copt418 * copt667;
  Real copt1225 = copt1221 + copt1222 + copt1223 + copt1224;
  Real copt1226 = -(copt1052 * copt1225 * copt681 * copt689);
  Real copt1227 = copt1220 + copt1226;
  Real copt940  = copt304 * copt438 * copt457;
  Real copt1240 = 2 * copt416 * copt457 * copt57;
  Real copt1241 = copt438 * copt459 * copt57;
  Real copt1242 = copt1240 + copt1241 + copt940;
  Real copt1257 = copt162 * copt449;
  Real copt1249 = copt232 + copt492;
  Real copt1250 = copt1249 * copt491;
  Real copt1251 = copt502 * copt682;
  Real copt1252 = copt443 + copt517;
  Real copt1253 = copt1252 * copt530;
  Real copt1254 = copt450 * copt538;
  Real copt1255 = copt1250 + copt1251 + copt1253 + copt1254;
  Real copt1256 = -(copt1016 * copt1255 * copt548 * copt577);
  Real copt1258 = -(copt449 * copt492);
  Real copt1259 = copt443 * copt558;
  Real copt1260 = -(copt162 * copt516);
  Real copt1261 = copt306 * copt516;
  Real copt1262 = copt1257 + copt1258 + copt1259 + copt1260 + copt1261;
  Real copt1263 = copt1262 * copt548;
  Real copt1264 = -(copt1022 * copt543 * copt577);
  Real copt1265 = copt1263 + copt1264;
  Real copt1266 = copt1016 * copt1265 * copt540;
  Real copt1267 = copt1256 + copt1266;
  Real copt1269 = copt588 * copt591;
  Real copt1270 = copt308 * copt598;
  Real copt1271 = copt604 * copt613;
  Real copt1272 = copt619 * copt686;
  Real copt1273 = copt1269 + copt1270 + copt1271 + copt1272;
  Real copt1274 = -(copt1035 * copt1273 * copt626 * copt647);
  Real copt1275 = copt443 * copt591;
  Real copt1276 = -(copt162 * copt602);
  Real copt1277 = copt1257 + copt1275 + copt1276 + copt600 + copt606;
  Real copt1278 = copt1277 * copt626;
  Real copt1279 = copt1211 * copt488 * copt647;
  Real copt1280 = copt1278 + copt1279;
  Real copt1281 = copt1035 * copt1280 * copt621;
  Real copt1282 = copt1274 + copt1281;
  Real copt1284 = copt308 * copt663;
  Real copt1285 = copt679 * copt686;
  Real copt1286 = copt1284 + copt1285;
  Real copt1287 = copt1052 * copt1286 * copt689 * copt715;
  Real copt1288 = -(copt162 * copt449);
  Real copt1289 = copt162 * copt667;
  Real copt1290 = copt1288 + copt1289 + copt665 + copt666 + copt668;
  Real copt1291 = -(copt1052 * copt1290 * copt681 * copt689);
  Real copt1292 = copt1287 + copt1291;
  Real copt950  = copt304 * copt450 * copt457;
  Real copt1305 = 2 * copt447 * copt457 * copt57;
  Real copt1306 = copt450 * copt459 * copt57;
  Real copt1307 = copt1305 + copt1306 + copt950;
  Real copt1322 = -(copt162 * copt418);
  Real copt1314 = copt162 + copt493;
  Real copt1315 = copt1314 * copt512;
  Real copt1316 = copt399 + copt497;
  Real copt1317 = copt1316 * copt530;
  Real copt1318 = copt308 * copt521;
  Real copt1319 = copt538 * copt684;
  Real copt1320 = copt1315 + copt1317 + copt1318 + copt1319;
  Real copt1321 = -(copt1016 * copt1320 * copt548 * copt577);
  Real copt1323 = copt317 * copt569;
  Real copt1324 = copt418 * copt492;
  Real copt1325 = copt162 * copt497;
  Real copt1326 = -(copt306 * copt497);
  Real copt1327 = copt1322 + copt1323 + copt1324 + copt1325 + copt1326;
  Real copt1328 = copt1327 * copt548;
  Real copt1329 = -(copt1022 * copt545 * copt577);
  Real copt1330 = copt1328 + copt1329;
  Real copt1331 = copt1016 * copt1330 * copt540;
  Real copt1332 = copt1321 + copt1331;
  Real copt1334 = copt512 * copt591;
  Real copt1335 = copt613 * copt615;
  Real copt1336 = copt308 * copt607;
  Real copt1337 = copt438 * copt619;
  Real copt1338 = copt1334 + copt1335 + copt1336 + copt1337;
  Real copt1339 = -(copt1035 * copt1338 * copt626 * copt647);
  Real copt1340 = copt317 * copt639;
  Real copt1341 = copt418 * copt589;
  Real copt1342 = copt162 * copt593;
  Real copt1343 = -(copt306 * copt593);
  Real copt1344 = copt1322 + copt1340 + copt1341 + copt1342 + copt1343;
  Real copt1345 = copt1344 * copt626;
  Real copt1346 = copt1211 * copt527 * copt647;
  Real copt1347 = copt1345 + copt1346;
  Real copt1348 = copt1035 * copt1347 * copt621;
  Real copt1349 = copt1339 + copt1348;
  Real copt1351 = copt308 * copt671;
  Real copt1352 = copt438 * copt679;
  Real copt1353 = copt1351 + copt1352;
  Real copt1354 = copt1052 * copt1353 * copt689 * copt715;
  Real copt1355 = copt162 * copt418;
  Real copt1356 = -(copt418 * copt654);
  Real copt1357 = copt317 * copt707;
  Real copt1358 = -(copt162 * copt659);
  Real copt1359 = copt306 * copt659;
  Real copt1360 = copt1355 + copt1356 + copt1357 + copt1358 + copt1359;
  Real copt1361 = -(copt1052 * copt1360 * copt681 * copt689);
  Real copt1362 = copt1354 + copt1361;
  Real copt960  = copt301 * copt304 * copt457;
  Real copt961  = 2 * copt304 * copt308;
  Real copt962  = copt302 + copt961;
  Real copt963  = copt459 * copt962;
  Real copt964  = copt960 + copt963;
  Real copt1393 = copt499 + copt593;
  Real copt1381 = -(copt446 * copt497);
  Real copt1382 = copt443 * copt500;
  Real copt1383 = copt317 * copt518;
  Real copt1384 = copt410 * copt516;
  Real copt1385 = copt1381 + copt1382 + copt1383 + copt1384;
  Real copt1386 = copt1016 * copt1385 * copt540 * copt548;
  Real copt1387 = copt502 * copt543;
  Real copt1388 = copt447 * copt521;
  Real copt1389 = copt1387 + copt1388;
  Real copt1390 = -(copt1016 * copt1389 * copt548 * copt577);
  Real copt1391 = copt1386 + copt1390;
  Real copt1394 = copt1393 * copt588;
  Real copt1395 = copt416 * copt598;
  Real copt1396 = copt509 + copt602;
  Real copt1397 = copt1396 * copt512;
  Real copt1398 = copt447 * copt607;
  Real copt1399 = copt1394 + copt1395 + copt1397 + copt1398;
  Real copt1400 = -(copt1035 * copt1399 * copt626 * copt647);
  Real copt1401 = -(copt446 * copt593);
  Real copt1402 = copt1393 * copt443;
  Real copt1403 = copt446 + copt603;
  Real copt1404 = copt1403 * copt317;
  Real copt1405 = copt410 * copt602;
  Real copt1406 = copt1401 + copt1402 + copt1404 + copt1405;
  Real copt1407 = copt1406 * copt626;
  Real copt1408 = -(copt1211 * copt505 * copt647);
  Real copt1409 = copt1407 + copt1408;
  Real copt1410 = copt1035 * copt1409 * copt621;
  Real copt1411 = copt1400 + copt1410;
  Real copt1413 = copt317 + copt673;
  Real copt1414 = copt1413 * copt588;
  Real copt1415 = copt416 * copt663;
  Real copt1416 = copt443 + copt709;
  Real copt1417 = copt1416 * copt512;
  Real copt1418 = copt447 * copt671;
  Real copt1419 = copt1414 + copt1415 + copt1417 + copt1418;
  Real copt1420 = copt1052 * copt1419 * copt689 * copt715;
  Real copt1421 = copt410 + copt673;
  Real copt1422 = copt1421 * copt443;
  Real copt1423 = copt446 * copt659;
  Real copt1424 = -(copt410 * copt667);
  Real copt1425 = copt509 + copt667;
  Real copt1426 = copt1425 * copt317;
  Real copt1427 = copt1422 + copt1423 + copt1424 + copt1426;
  Real copt1428 = -(copt1427 * copt689);
  Real copt1429 = copt1058 * copt682 * copt715;
  Real copt1430 = copt1428 + copt1429;
  Real copt1431 = copt1052 * copt1430 * copt681;
  Real copt1432 = copt1420 + copt1431;
  Real copt970  = copt440 * copt459;
  Real copt971  = copt304 * copt465;
  Real copt972  = copt970 + copt971;
  Real copt1451 = -(copt162 * copt446);
  Real copt1464 = copt300 + copt638;
  Real copt1452 = copt443 * copt494;
  Real copt1453 = copt446 * copt492;
  Real copt1454 = copt162 * copt516;
  Real copt1455 = -(copt300 * copt516);
  Real copt1456 = copt1451 + copt1452 + copt1453 + copt1454 + copt1455;
  Real copt1457 = copt1016 * copt1456 * copt540 * copt548;
  Real copt1458 = copt301 * copt502;
  Real copt1459 = copt538 * copt545;
  Real copt1460 = copt1458 + copt1459;
  Real copt1461 = -(copt1016 * copt1460 * copt548 * copt577);
  Real copt1462 = copt1457 + copt1461;
  Real copt1465 = copt1464 * copt588;
  Real copt1466 = copt541 * copt598;
  Real copt1467 = copt1396 * copt613;
  Real copt1468 = copt447 * copt619;
  Real copt1469 = copt1465 + copt1466 + copt1467 + copt1468;
  Real copt1470 = -(copt1035 * copt1469 * copt626 * copt647);
  Real copt1471 = copt1464 * copt443;
  Real copt1472 = copt446 * copt589;
  Real copt1473 = copt162 * copt602;
  Real copt1474 = -(copt300 * copt602);
  Real copt1475 = copt1451 + copt1471 + copt1472 + copt1473 + copt1474;
  Real copt1476 = copt1475 * copt626;
  Real copt1477 = -(copt1211 * copt488 * copt647);
  Real copt1478 = copt1476 + copt1477;
  Real copt1479 = copt1035 * copt1478 * copt621;
  Real copt1480 = copt1470 + copt1479;
  Real copt1482 = copt232 + copt654;
  Real copt1483 = copt1482 * copt588;
  Real copt1484 = copt541 * copt663;
  Real copt1485 = copt1416 * copt613;
  Real copt1486 = copt447 * copt679;
  Real copt1487 = copt1483 + copt1484 + copt1485 + copt1486;
  Real copt1488 = copt1052 * copt1487 * copt689 * copt715;
  Real copt1489 = copt162 * copt446;
  Real copt1490 = -(copt446 * copt654);
  Real copt1491 = copt484 + copt654;
  Real copt1492 = copt1491 * copt443;
  Real copt1493 = -(copt162 * copt667);
  Real copt1494 = copt300 * copt667;
  Real copt1495 = copt1489 + copt1490 + copt1492 + copt1493 + copt1494;
  Real copt1496 = -(copt1495 * copt689);
  Real copt1497 = copt1058 * copt684 * copt715;
  Real copt1498 = copt1496 + copt1497;
  Real copt1499 = copt1052 * copt1498 * copt681;
  Real copt1500 = copt1488 + copt1499;
  Real copt977  = copt452 * copt459;
  Real copt978  = copt304 * copt469;
  Real copt979  = copt977 + copt978;
  Real copt1519 = copt162 * copt410;
  Real copt1520 = -(copt410 * copt492);
  Real copt1521 = copt317 * copt514;
  Real copt1522 = -(copt162 * copt497);
  Real copt1523 = copt300 * copt497;
  Real copt1524 = copt1519 + copt1520 + copt1521 + copt1522 + copt1523;
  Real copt1525 = copt1016 * copt1524 * copt540 * copt548;
  Real copt1526 = copt521 * copt541;
  Real copt1527 = copt416 * copt538;
  Real copt1528 = copt1526 + copt1527;
  Real copt1529 = -(copt1016 * copt1528 * copt548 * copt577);
  Real copt1530 = copt1525 + copt1529;
  Real copt1532 = copt1464 * copt512;
  Real copt1533 = copt410 + copt594;
  Real copt1534 = copt1533 * copt613;
  Real copt1535 = copt541 * copt607;
  Real copt1536 = copt543 * copt619;
  Real copt1537 = copt1532 + copt1534 + copt1535 + copt1536;
  Real copt1538 = -(copt1035 * copt1537 * copt626 * copt647);
  Real copt1539 = -(copt410 * copt589);
  Real copt1540 = copt484 + copt589;
  Real copt1541 = copt1540 * copt317;
  Real copt1542 = -(copt162 * copt593);
  Real copt1543 = copt300 * copt593;
  Real copt1544 = copt1519 + copt1539 + copt1541 + copt1542 + copt1543;
  Real copt1545 = copt1544 * copt626;
  Real copt1546 = -(copt1211 * copt527 * copt647);
  Real copt1547 = copt1545 + copt1546;
  Real copt1548 = copt1035 * copt1547 * copt621;
  Real copt1549 = copt1538 + copt1548;
  Real copt1551 = copt1482 * copt512;
  Real copt1552 = copt399 + copt659;
  Real copt1553 = copt1552 * copt613;
  Real copt1554 = copt541 * copt671;
  Real copt1555 = copt543 * copt679;
  Real copt1556 = copt1551 + copt1553 + copt1554 + copt1555;
  Real copt1557 = copt1052 * copt1556 * copt689 * copt715;
  Real copt1558 = -(copt162 * copt410);
  Real copt1559 = copt300 + copt655;
  Real copt1560 = copt1559 * copt317;
  Real copt1561 = copt410 * copt654;
  Real copt1562 = copt162 * copt659;
  Real copt1563 = -(copt300 * copt659);
  Real copt1564 = copt1558 + copt1560 + copt1561 + copt1562 + copt1563;
  Real copt1565 = -(copt1564 * copt689);
  Real copt1566 = copt1058 * copt686 * copt715;
  Real copt1567 = copt1565 + copt1566;
  Real copt1568 = copt1052 * copt1567 * copt681;
  Real copt1569 = copt1557 + copt1568;
  Real copt1581 = copt1016 * copt540 * copt548 * copt613;
  Real copt1582 = copt416 * copt491;
  Real copt1583 = copt512 * copt545;
  Real copt1584 = copt1582 + copt1583;
  Real copt1585 = -(copt1016 * copt1584 * copt548 * copt577);
  Real copt1586 = copt1581 + copt1585;
  Real copt1591 = -(copt306 * copt446);
  Real copt1592 = copt443 * copt485;
  Real copt1593 = copt300 * copt449;
  Real copt1594 = copt1288 + copt1489 + copt1591 + copt1592 + copt1593;
  Real copt1595 = copt1016 * copt1594 * copt540 * copt548;
  Real copt1596 = copt491 * copt541;
  Real copt1597 = copt447 * copt530;
  Real copt1598 = copt1596 + copt1597;
  Real copt1599 = -(copt1016 * copt1598 * copt548 * copt577);
  Real copt1600 = copt1595 + copt1599;
  Real copt1605 = copt1355 + copt1558 + copt584 + copt585 + copt586;
  Real copt1606 = copt1016 * copt1605 * copt540 * copt548;
  Real copt1607 = copt530 * copt543;
  Real copt1608 = copt301 * copt512;
  Real copt1609 = copt1607 + copt1608;
  Real copt1610 = -(copt1016 * copt1609 * copt548 * copt577);
  Real copt1611 = copt1606 + copt1610;
  Real copt1616 = copt1035 * copt613 * copt621 * copt626;
  Real copt1617 = copt488 * copt588;
  Real copt1618 = copt512 * copt527;
  Real copt1619 = copt1617 + copt1618;
  Real copt1620 = -(copt1035 * copt1619 * copt626 * copt647);
  Real copt1621 = copt1616 + copt1620;
  Real copt1626 = copt1035 * copt1594 * copt621 * copt626;
  Real copt1627 = copt485 * copt588;
  Real copt1628 = copt527 * copt613;
  Real copt1629 = copt1627 + copt1628;
  Real copt1630 = -(copt1035 * copt1629 * copt626 * copt647);
  Real copt1631 = copt1626 + copt1630;
  Real copt1636 = copt1035 * copt1605 * copt621 * copt626;
  Real copt1637 = copt485 * copt512;
  Real copt1638 = copt524 * copt613;
  Real copt1639 = copt1637 + copt1638;
  Real copt1640 = -(copt1035 * copt1639 * copt626 * copt647);
  Real copt1641 = copt1636 + copt1640;
  Real copt1646 = copt438 * copt588;
  Real copt1647 = copt450 * copt512;
  Real copt1648 = copt1646 + copt1647;
  Real copt1649 = copt1052 * copt1648 * copt689 * copt715;
  Real copt1650 = -(copt1052 * copt530 * copt681 * copt689);
  Real copt1651 = copt1649 + copt1650;
  Real copt1656 = copt588 * copt682;
  Real copt1657 = copt450 * copt613;
  Real copt1658 = copt1656 + copt1657;
  Real copt1659 = copt1052 * copt1658 * copt689 * copt715;
  Real copt1660 = copt1257 + copt1451 + copt506 + copt507 + copt508;
  Real copt1661 = -(copt1052 * copt1660 * copt681 * copt689);
  Real copt1662 = copt1659 + copt1661;
  Real copt1667 = copt512 * copt682;
  Real copt1668 = copt613 * copt684;
  Real copt1669 = copt1667 + copt1668;
  Real copt1670 = copt1052 * copt1669 * copt689 * copt715;
  Real copt1671 = copt1322 + copt1519 + copt483 + copt486 + copt490;
  Real copt1672 = -(copt1052 * copt1671 * copt681 * copt689);
  Real copt1673 = copt1670 + copt1672;
  Real copt1695 = Power(copt806, 2);
  Real copt1696 = 1 / copt1695;
  Real copt820  = copt815 * copt819;
  Real copt831  = 2 * copt457 * copt459 * copt830;
  Real copt837  = copt832 * copt836;
  Real copt838  = copt820 + copt831 + copt837;
  Real copt839  = -(copt814 * copt838);
  Real copt859  = copt849 * copt858;
  Real copt860  = copt839 + copt859;
  Real copt1698 = 2 * copt162;
  Real copt1702 = copt1698 + copt484 + copt504;
  Real copt1699 = -2 * copt306;
  Real copt1700 = copt1698 + copt1699;
  Real copt1066 = copt1026 * copt581 * copt582 * l1 * l2;
  Real copt1067 = copt1041 * copt651 * copt652 * l0 * l2;
  Real copt1068 = copt1062 * copt719 * copt720 * l0 * l1;
  Real copt1069 = copt1066 + copt1067 + copt1068;
  Real copt1732 = 2 * copt317;
  Real copt1736 = copt1732 + copt487 + copt499;
  Real copt1733 = -2 * copt418;
  Real copt1734 = copt1732 + copt1733;
  Real copt1114 = copt1091 * copt581 * copt582 * l1 * l2;
  Real copt1115 = copt1098 * copt651 * copt652 * l0 * l2;
  Real copt1116 = copt1110 * copt719 * copt720 * l0 * l1;
  Real copt1117 = copt1114 + copt1115 + copt1116;
  Real copt1767 = 2 * copt443;
  Real copt1771 = copt1767 + copt509 + copt526;
  Real copt1768 = -2 * copt449;
  Real copt1769 = copt1767 + copt1768;
  Real copt1162 = copt1139 * copt581 * copt582 * l1 * l2;
  Real copt1163 = copt1146 * copt651 * copt652 * l0 * l2;
  Real copt1164 = copt1158 * copt719 * copt720 * l0 * l1;
  Real copt1165 = copt1162 + copt1163 + copt1164;
  Real copt1797 = 2 * copt300;
  Real copt1798 = copt1699 + copt1797;
  Real copt1820 = -2 * copt162;
  Real copt1821 = copt1797 + copt1820;
  Real copt1231 = copt1197 * copt581 * copt582 * l1 * l2;
  Real copt1232 = copt1215 * copt651 * copt652 * l0 * l2;
  Real copt1233 = copt1227 * copt719 * copt720 * l0 * l1;
  Real copt1234 = copt1231 + copt1232 + copt1233;
  Real copt1843 = 2 * copt410;
  Real copt1864 = -2 * copt317;
  Real copt1865 = copt1843 + copt1864;
  Real copt1296 = copt1267 * copt581 * copt582 * l1 * l2;
  Real copt1297 = copt1282 * copt651 * copt652 * l0 * l2;
  Real copt1298 = copt1292 * copt719 * copt720 * l0 * l1;
  Real copt1299 = copt1296 + copt1297 + copt1298;
  Real copt1856 = copt162 * copt505;
  Real copt1857 = -(copt300 * copt306);
  Real copt1853 = -(copt418 * copt449);
  Real copt1891 = 2 * copt446;
  Real copt1904 = -2 * copt443;
  Real copt1905 = copt1891 + copt1904;
  Real copt1366 = copt1332 * copt581 * copt582 * l1 * l2;
  Real copt1367 = copt1349 * copt651 * copt652 * l0 * l2;
  Real copt1368 = copt1362 * copt719 * copt720 * l0 * l1;
  Real copt1369 = copt1366 + copt1367 + copt1368;
  Real copt1712 = -2 * copt300;
  Real copt1928 = 2 * copt306;
  Real copt1929 = copt1712 + copt1928;
  Real copt1801 = 2 * copt162 * copt410 * copt418;
  Real copt1809 = 2 * copt162 * copt446 * copt449;
  Real copt1947 = copt1820 + copt1928;
  Real copt1436 = copt1391 * copt581 * copt582 * l1 * l2;
  Real copt1437 = copt1411 * copt651 * copt652 * l0 * l2;
  Real copt1438 = copt1432 * copt719 * copt720 * l0 * l1;
  Real copt1439 = copt1436 + copt1437 + copt1438;
  Real copt1746 = -2 * copt410;
  Real copt1899 = 2 * copt418 * copt446;
  Real copt1859 = -(copt446 * copt449);
  Real copt1974 = 2 * copt418;
  Real copt1988 = copt1864 + copt1974;
  Real copt1504 = copt1462 * copt581 * copt582 * l1 * l2;
  Real copt1505 = copt1480 * copt651 * copt652 * l0 * l2;
  Real copt1506 = copt1500 * copt719 * copt720 * l0 * l1;
  Real copt1507 = copt1504 + copt1505 + copt1506;
  Real copt1981 = -(copt162 * copt505);
  Real copt1888 = -(copt410 * copt418);
  Real copt1781 = -2 * copt446;
  Real copt1978 = -(copt410 * copt446);
  Real copt1852 = 2 * copt410 * copt449;
  Real copt2018 = 2 * copt449;
  Real copt2027 = copt1904 + copt2018;
  Real copt1573 = copt1530 * copt581 * copt582 * l1 * l2;
  Real copt1574 = copt1549 * copt651 * copt652 * l0 * l2;
  Real copt1575 = copt1569 * copt719 * copt720 * l0 * l1;
  Real copt1576 = copt1573 + copt1574 + copt1575;
  Real copt1678 = 2 * copt162 * copt745;
  Real copt1679 = 2 * copt162 * copt747;
  Real copt1680 = -2 * copt306 * copt745;
  Real copt1681 = -2 * copt306 * copt747;
  Real copt1682 = -2 * copt317 * copt488 * copt505;
  Real copt1683 = -4 * copt162 * copt410 * copt418;
  Real copt1684 = 2 * copt300 * copt410 * copt418;
  Real copt1685 = 2 * copt306 * copt410 * copt418;
  Real copt1686 = 2 * copt162 * copt758;
  Real copt1687 = -2 * copt300 * copt758;
  Real copt1688 = -2 * copt443 * copt505 * copt527;
  Real copt1689 = -4 * copt162 * copt446 * copt449;
  Real copt1690 = 2 * copt300 * copt446 * copt449;
  Real copt1691 = 2 * copt306 * copt446 * copt449;
  Real copt1692 = 2 * copt162 * copt774;
  Real copt1693 = -2 * copt300 * copt774;
  Real copt1694 = copt1678 + copt1679 + copt1680 + copt1681 + copt1682 +
                  copt1683 + copt1684 + copt1685 + copt1686 + copt1687 +
                  copt1688 + copt1689 + copt1690 + copt1691 + copt1692 +
                  copt1693;
  Real copt863  = copt819 * copt862;
  Real copt864  = 2 * copt304 * copt57 * copt830;
  Real copt866  = copt836 * copt865;
  Real copt867  = copt863 + copt864 + copt866;
  Real copt868  = -(copt849 * copt867);
  Real copt869  = copt814 * copt858;
  Real copt870  = copt868 + copt869;
  Real copt1701 = copt1700 * copt459;
  Real copt1703 = copt1702 * copt457;
  Real copt1704 = copt1701 + copt1703;
  Real copt1705 = copt1704 * copt304;
  Real copt1706 = 2 * copt457 * copt541;
  Real copt1707 = copt1702 * copt459;
  Real copt1708 = copt1706 + copt1707;
  Real copt1709 = copt1708 * copt57;
  Real copt1710 = copt1705 + copt1709;
  Real copt1713 = copt1698 + copt1712;
  Real copt1720 = copt1026 * copt840 * l1 * l2;
  Real copt1721 = copt1041 * copt843 * l0 * l2;
  Real copt1722 = copt1062 * copt846 * l0 * l1;
  Real copt1723 = copt1720 + copt1721 + copt1722;
  Real copt1728 = 2 * copt317 * copt781;
  Real copt1729 = -2 * copt793;
  Real copt1730 = copt1728 + copt1729;
  Real copt1735 = copt1734 * copt459;
  Real copt1737 = copt1736 * copt457;
  Real copt1738 = copt1735 + copt1737;
  Real copt1739 = copt1738 * copt304;
  Real copt1740 = 2 * copt457 * copt543;
  Real copt1741 = copt1736 * copt459;
  Real copt1742 = copt1740 + copt1741;
  Real copt1743 = copt1742 * copt57;
  Real copt1744 = copt1739 + copt1743;
  Real copt1747 = copt1732 + copt1746;
  Real copt1754 = copt1091 * copt840 * l1 * l2;
  Real copt1755 = copt1098 * copt843 * l0 * l2;
  Real copt1756 = copt1110 * copt846 * l0 * l1;
  Real copt1757 = copt1754 + copt1755 + copt1756;
  Real copt1762 = 2 * copt443 * copt767;
  Real copt1763 = -2 * copt317 * copt488 * copt527;
  Real copt1764 = -2 * copt804;
  Real copt1765 = copt1762 + copt1763 + copt1764;
  Real copt1770 = copt1769 * copt459;
  Real copt1772 = copt1771 * copt457;
  Real copt1773 = copt1770 + copt1772;
  Real copt1774 = copt1773 * copt304;
  Real copt1775 = 2 * copt457 * copt545;
  Real copt1776 = copt1771 * copt459;
  Real copt1777 = copt1775 + copt1776;
  Real copt1778 = copt1777 * copt57;
  Real copt1779 = copt1774 + copt1778;
  Real copt1782 = copt1767 + copt1781;
  Real copt1789 = copt1139 * copt840 * l1 * l2;
  Real copt1790 = copt1146 * copt843 * l0 * l2;
  Real copt1791 = copt1158 * copt846 * l0 * l1;
  Real copt1792 = copt1789 + copt1790 + copt1791;
  Real copt1799 = copt1798 * copt779;
  Real copt1800 = copt1798 * copt764;
  Real copt1802 = -2 * copt306 * copt410 * copt418;
  Real copt1803 = -2 * copt162 * copt758;
  Real copt1804 = 2 * copt300 * copt758;
  Real copt1805 = 2 * copt300 * copt418;
  Real copt1806 = -(copt306 * copt418);
  Real copt1807 = copt1805 + copt1806 + copt483 + copt489;
  Real copt1808 = -2 * copt1807 * copt317;
  Real copt1810 = -2 * copt306 * copt446 * copt449;
  Real copt1811 = -2 * copt162 * copt774;
  Real copt1812 = 2 * copt300 * copt774;
  Real copt1813 = copt162 * copt527;
  Real copt1814 = 2 * copt300 * copt449;
  Real copt1815 = -(copt306 * copt802);
  Real copt1816 = copt1813 + copt1814 + copt1815;
  Real copt1817 = -2 * copt1816 * copt443;
  Real copt1818 = copt1799 + copt1800 + copt1801 + copt1802 + copt1803 +
                  copt1804 + copt1808 + copt1809 + copt1810 + copt1811 +
                  copt1812 + copt1817;
  Real copt1826 = copt1821 * copt457;
  Real copt1827 = copt1826 + copt460;
  Real copt1828 = copt1827 * copt57;
  Real copt1829 = copt1828 + copt930;
  Real copt1832 = copt1197 * copt840 * l1 * l2;
  Real copt1833 = copt1215 * copt843 * l0 * l2;
  Real copt1834 = copt1227 * copt846 * l0 * l1;
  Real copt1835 = copt1832 + copt1833 + copt1834;
  Real copt1840 = 2 * copt410 * copt744;
  Real copt1841 = -4 * copt162 * copt306 * copt410;
  Real copt1842 = 2 * copt410 * copt751;
  Real copt1844 = copt1733 + copt1843;
  Real copt1845 = copt1844 * copt764;
  Real copt1846 = -2 * copt418 * copt744;
  Real copt1847 = 2 * copt162 * copt300 * copt418;
  Real copt1848 = 2 * copt162 * copt306 * copt418;
  Real copt1849 = -2 * copt300 * copt306 * copt418;
  Real copt1850 = -2 * copt418 * copt446 * copt449;
  Real copt1851 = 2 * copt410 * copt774;
  Real copt1854 = copt1852 + copt1853 + copt523;
  Real copt1855 = -2 * copt1854 * copt443;
  Real copt1858 = copt443 * copt527;
  Real copt1860 = copt1856 + copt1857 + copt1858 + copt1859 + copt751 + copt774;
  Real copt1861 = -2 * copt1860 * copt317;
  Real copt1862 = copt1840 + copt1841 + copt1842 + copt1845 + copt1846 +
                  copt1847 + copt1848 + copt1849 + copt1850 + copt1851 +
                  copt1855 + copt1861;
  Real copt1870 = copt1865 * copt457;
  Real copt1871 = copt1870 + copt464;
  Real copt1872 = copt1871 * copt57;
  Real copt1873 = copt1872 + copt940;
  Real copt1876 = copt1267 * copt840 * l1 * l2;
  Real copt1877 = copt1282 * copt843 * l0 * l2;
  Real copt1878 = copt1292 * copt846 * l0 * l1;
  Real copt1879 = copt1876 + copt1877 + copt1878;
  Real copt1884 = 2 * copt446 * copt744;
  Real copt1885 = -4 * copt162 * copt306 * copt446;
  Real copt1886 = 2 * copt446 * copt751;
  Real copt1887 = 2 * copt446 * copt758;
  Real copt1889 = copt1856 + copt1857 + copt1888 + copt751 + copt758;
  Real copt1890 = -2 * copt1889 * copt443;
  Real copt1892 = copt1768 + copt1891;
  Real copt1893 = copt1892 * copt779;
  Real copt1894 = -2 * copt449 * copt744;
  Real copt1895 = 2 * copt162 * copt300 * copt449;
  Real copt1896 = 2 * copt162 * copt306 * copt449;
  Real copt1897 = -2 * copt300 * copt306 * copt449;
  Real copt1898 = -2 * copt410 * copt418 * copt449;
  Real copt1900 = copt1853 + copt1899 + copt609 + copt611;
  Real copt1901 = -2 * copt1900 * copt317;
  Real copt1902 = copt1884 + copt1885 + copt1886 + copt1887 + copt1890 +
                  copt1893 + copt1894 + copt1895 + copt1896 + copt1897 +
                  copt1898 + copt1901;
  Real copt1910 = copt1905 * copt457;
  Real copt1911 = copt1910 + copt468;
  Real copt1912 = copt1911 * copt57;
  Real copt1913 = copt1912 + copt950;
  Real copt1916 = copt1332 * copt840 * l1 * l2;
  Real copt1917 = copt1349 * copt843 * l0 * l2;
  Real copt1918 = copt1362 * copt846 * l0 * l1;
  Real copt1919 = copt1916 + copt1917 + copt1918;
  Real copt1924 = -2 * copt162 * copt745;
  Real copt1925 = -2 * copt162 * copt747;
  Real copt1926 = 2 * copt306 * copt745;
  Real copt1927 = 2 * copt306 * copt747;
  Real copt1930 = copt1929 * copt779;
  Real copt1931 = copt1929 * copt764;
  Real copt1932 = -2 * copt300 * copt410 * copt418;
  Real copt1933 = -(copt300 * copt410);
  Real copt1934 = 2 * copt306 * copt410;
  Real copt1935 = -(copt162 * copt488);
  Real copt1936 = copt1933 + copt1934 + copt1935 + copt586;
  Real copt1937 = -2 * copt1936 * copt317;
  Real copt1938 = -2 * copt300 * copt446 * copt449;
  Real copt1939 = 2 * copt306 * copt446;
  Real copt1940 = -(copt162 * copt527);
  Real copt1941 = -(copt300 * copt802);
  Real copt1942 = copt1939 + copt1940 + copt1941;
  Real copt1943 = -2 * copt1942 * copt443;
  Real copt1944 = copt1801 + copt1809 + copt1924 + copt1925 + copt1926 +
                  copt1927 + copt1930 + copt1931 + copt1932 + copt1937 +
                  copt1938 + copt1943;
  Real copt1951 = copt301 * copt459 * copt57;
  Real copt1952 = copt1947 * copt459;
  Real copt1953 = copt1952 + copt458;
  Real copt1954 = copt1953 * copt304;
  Real copt1955 = copt1951 + copt1954;
  Real copt1958 = copt1391 * copt840 * l1 * l2;
  Real copt1959 = copt1411 * copt843 * l0 * l2;
  Real copt1960 = copt1432 * copt846 * l0 * l1;
  Real copt1961 = copt1958 + copt1959 + copt1960;
  Real copt1966 = -2 * copt410 * copt744;
  Real copt1967 = 2 * copt162 * copt300 * copt410;
  Real copt1968 = 2 * copt162 * copt306 * copt410;
  Real copt1969 = -2 * copt300 * copt306 * copt410;
  Real copt1970 = 2 * copt418 * copt744;
  Real copt1971 = -4 * copt162 * copt300 * copt418;
  Real copt1972 = 2 * copt418 * copt761;
  Real copt1973 = 2 * copt418 * copt747;
  Real copt1975 = copt1746 + copt1974;
  Real copt1976 = copt1975 * copt764;
  Real copt1977 = -2 * copt410 * copt446 * copt449;
  Real copt1979 = copt1899 + copt1978 + copt611;
  Real copt1980 = -2 * copt1979 * copt443;
  Real copt1982 = -(copt443 * copt527);
  Real copt1983 = copt1857 + copt1859 + copt1981 + copt1982 + copt747 + copt761;
  Real copt1984 = -2 * copt1983 * copt317;
  Real copt1985 = copt1966 + copt1967 + copt1968 + copt1969 + copt1970 +
                  copt1971 + copt1972 + copt1973 + copt1976 + copt1977 +
                  copt1980 + copt1984;
  Real copt1992 = copt416 * copt459 * copt57;
  Real copt1993 = copt1988 * copt459;
  Real copt1994 = copt1993 + copt463;
  Real copt1995 = copt1994 * copt304;
  Real copt1996 = copt1992 + copt1995;
  Real copt1999 = copt1462 * copt840 * l1 * l2;
  Real copt2000 = copt1480 * copt843 * l0 * l2;
  Real copt2001 = copt1500 * copt846 * l0 * l1;
  Real copt2002 = copt1999 + copt2000 + copt2001;
  Real copt2007 = -2 * copt446 * copt744;
  Real copt2008 = 2 * copt162 * copt300 * copt446;
  Real copt2009 = 2 * copt162 * copt306 * copt446;
  Real copt2010 = -2 * copt300 * copt306 * copt446;
  Real copt2011 = -2 * copt410 * copt418 * copt446;
  Real copt2012 = copt1857 + copt1888 + copt1981 + copt745 + copt761;
  Real copt2013 = -2 * copt2012 * copt443;
  Real copt2014 = 2 * copt449 * copt744;
  Real copt2015 = -4 * copt162 * copt300 * copt449;
  Real copt2016 = 2 * copt449 * copt761;
  Real copt2017 = 2 * copt449 * copt745;
  Real copt2019 = copt1781 + copt2018;
  Real copt2020 = copt2019 * copt779;
  Real copt2021 = -(copt443 * copt488);
  Real copt2022 = copt1852 + copt1978 + copt2021 + copt523;
  Real copt2023 = -2 * copt2022 * copt317;
  Real copt2024 = copt2007 + copt2008 + copt2009 + copt2010 + copt2011 +
                  copt2013 + copt2014 + copt2015 + copt2016 + copt2017 +
                  copt2020 + copt2023;
  Real copt2031 = copt447 * copt459 * copt57;
  Real copt2032 = copt2027 * copt459;
  Real copt2033 = copt2032 + copt467;
  Real copt2034 = copt2033 * copt304;
  Real copt2035 = copt2031 + copt2034;
  Real copt2038 = copt1530 * copt840 * l1 * l2;
  Real copt2039 = copt1549 * copt843 * l0 * l2;
  Real copt2040 = copt1569 * copt846 * l0 * l1;
  Real copt2041 = copt2038 + copt2039 + copt2040;
  out1(0)       = copt455;
  out1(1)       = copt456 * copt471 * copt477;
  out1(2)       = copt476;
  out1(3) = (copt479 * copt480 * copt481 * copt482 * copt732 * copt737) / 2.;
  out1(4) =
      (copt479 * copt480 * copt481 * copt482 * copt743 * copt807 * copt860) /
      2.;
  out1(5) =
      (copt479 * copt480 * copt481 * copt482 * copt743 * copt807 * copt870) /
      2.;
  out2(0, 0)  = copt456 * copt872 * copt875;
  out2(0, 1)  = copt456 * copt872 * copt879;
  out2(0, 2)  = copt456 * copt872 * copt883;
  out2(0, 3)  = copt310 * copt456 * copt57;
  out2(0, 4)  = copt440 * copt456 * copt57;
  out2(0, 5)  = copt452 * copt456 * copt57;
  out2(0, 6)  = copt304 * copt310 * copt456;
  out2(0, 7)  = copt304 * copt440 * copt456;
  out2(0, 8)  = copt304 * copt452 * copt456;
  out2(0, 9)  = 0;
  out2(0, 10) = 0;
  out2(0, 11) = 0;
  out2(0, 12) = 0;
  out2(0, 13) = 0;
  out2(0, 14) = 0;
  out2(0, 15) = 0;
  out2(0, 16) = 0;
  out2(0, 17) = 0;
  out2(1, 0) =
      copt892 * copt894 *
      (copt310 * copt471 * copt475 * copt872 +
       copt454 * copt461 * copt471 * copt895 + copt454 * copt475 * copt902);
  out2(1, 1) =
      copt892 * copt894 *
      (copt440 * copt471 * copt475 * copt872 +
       copt454 * copt465 * copt471 * copt895 + copt454 * copt475 * copt913);
  out2(1, 2) =
      copt892 * copt894 *
      (copt452 * copt471 * copt475 * copt872 +
       copt454 * copt469 * copt471 * copt895 + copt454 * copt475 * copt924);
  out2(1, 3) = copt892 * copt894 *
               (-(copt454 * copt457 * copt461 * copt471) -
                copt310 * copt471 * copt475 * copt57 +
                copt454 * copt475 *
                    ((copt460 - 2 * copt457 * copt541) * copt57 + copt930));
  out2(1, 4) = copt892 * copt894 *
               (-(copt454 * copt457 * copt465 * copt471) -
                copt440 * copt471 * copt475 * copt57 +
                copt454 * copt475 *
                    ((copt464 - 2 * copt457 * copt543) * copt57 + copt940));
  out2(1, 5) = copt892 * copt894 *
               (-(copt454 * copt457 * copt469 * copt471) -
                copt452 * copt471 * copt475 * copt57 +
                copt454 * copt475 *
                    ((copt468 - 2 * copt457 * copt545) * copt57 + copt950));
  out2(1, 6) =
      copt892 * copt894 *
      (-(copt454 * copt459 * copt461 * copt471) -
       copt304 * copt310 * copt471 * copt475 + copt454 * copt475 * copt964);
  out2(1, 7) = -(copt304 * copt440 * copt471 * copt477 * copt892) -
               copt456 * copt459 * copt465 * copt471 * copt894 +
               copt456 * copt477 * copt972;
  out2(1, 8) = -(copt304 * copt452 * copt471 * copt477 * copt892) -
               copt456 * copt459 * copt469 * copt471 * copt894 +
               copt456 * copt477 * copt979;
  out2(1, 9)  = 0;
  out2(1, 10) = 0;
  out2(1, 11) = 0;
  out2(1, 12) = 0;
  out2(1, 13) = 0;
  out2(1, 14) = 0;
  out2(1, 15) = 0;
  out2(1, 16) = 0;
  out2(1, 17) = 0;
  out2(2, 0)  = copt461 * copt477 * copt985;
  out2(2, 1)  = copt465 * copt477 * copt985;
  out2(2, 2)  = copt469 * copt477 * copt985;
  out2(2, 3)  = copt457 * copt461 * copt477;
  out2(2, 4)  = copt457 * copt465 * copt477;
  out2(2, 5)  = copt457 * copt469 * copt477;
  out2(2, 6)  = copt459 * copt461 * copt477;
  out2(2, 7)  = copt459 * copt465 * copt477;
  out2(2, 8)  = copt459 * copt469 * copt477;
  out2(2, 9)  = 0;
  out2(2, 10) = 0;
  out2(2, 11) = 0;
  out2(2, 12) = 0;
  out2(2, 13) = 0;
  out2(2, 14) = 0;
  out2(2, 15) = 0;
  out2(2, 16) = 0;
  out2(2, 17) = 0;
  out2(3, 0) =
      -(copt1003 * copt479 * copt480 * copt481 * copt482 * copt732 *
        (-2 * copt471 * copt902 + 2 * copt454 * copt461 * copt985 +
         2 * copt310 * copt475 * copt999)) /
          2. +
      (copt479 * copt480 * copt481 * copt482 * copt737 *
       (copt1069 * copt471 + copt722 * copt902 -
        2 * copt461 * copt730 * copt985 -
        copt475 * (copt1062 * copt728 * l0 * l1 + copt1041 * copt726 * l0 * l2 +
                   copt1026 * copt724 * l1 * l2))) /
          2.;
  out2(3, 1) =
      -(copt1003 * copt479 * copt480 * copt481 * copt482 * copt732 *
        (-2 * copt471 * copt913 + 2 * copt454 * copt465 * copt985 +
         2 * copt440 * copt475 * copt999)) /
          2. +
      (copt479 * copt480 * copt481 * copt482 * copt737 *
       (copt1117 * copt471 + copt722 * copt913 -
        2 * copt465 * copt730 * copt985 -
        copt475 * (copt1110 * copt728 * l0 * l1 + copt1098 * copt726 * l0 * l2 +
                   copt1091 * copt724 * l1 * l2))) /
          2.;
  out2(3, 2) =
      -(copt1003 * copt479 * copt480 * copt481 * copt482 * copt732 *
        (-2 * copt471 * copt924 + 2 * copt454 * copt469 * copt985 +
         2 * copt452 * copt475 * copt999)) /
          2. +
      (copt479 * copt480 * copt481 * copt482 * copt737 *
       (copt1165 * copt471 + copt722 * copt924 -
        2 * copt469 * copt730 * copt985 -
        copt475 * (copt1158 * copt728 * l0 * l1 + copt1146 * copt726 * l0 * l2 +
                   copt1139 * copt724 * l1 * l2))) /
          2.;
  out2(3, 3) =
      -(copt1003 * copt479 * copt480 * copt481 * copt482 *
        (2 * copt454 * copt457 * copt461 - 2 * copt1173 * copt471 +
         2 * copt310 * copt475 * copt57) *
        copt732) /
          2. +
      (copt479 * copt480 * copt481 * copt482 * copt737 *
       (copt1234 * copt471 + copt1173 * copt722 -
        2 * copt457 * copt461 * copt730 -
        copt475 * (copt1227 * copt728 * l0 * l1 + copt1215 * copt726 * l0 * l2 +
                   copt1197 * copt724 * l1 * l2))) /
          2.;
  out2(3, 4) =
      -(copt1003 * copt479 * copt480 * copt481 * copt482 *
        (2 * copt454 * copt457 * copt465 - 2 * copt1242 * copt471 +
         2 * copt440 * copt475 * copt57) *
        copt732) /
          2. +
      (copt479 * copt480 * copt481 * copt482 * copt737 *
       (copt1299 * copt471 + copt1242 * copt722 -
        2 * copt457 * copt465 * copt730 -
        copt475 * (copt1292 * copt728 * l0 * l1 + copt1282 * copt726 * l0 * l2 +
                   copt1267 * copt724 * l1 * l2))) /
          2.;
  out2(3, 5) =
      -(copt1003 * copt479 * copt480 * copt481 * copt482 *
        (2 * copt454 * copt457 * copt469 - 2 * copt1307 * copt471 +
         2 * copt452 * copt475 * copt57) *
        copt732) /
          2. +
      (copt479 * copt480 * copt481 * copt482 * copt737 *
       (copt1369 * copt471 + copt1307 * copt722 -
        2 * copt457 * copt469 * copt730 -
        copt475 * (copt1362 * copt728 * l0 * l1 + copt1349 * copt726 * l0 * l2 +
                   copt1332 * copt724 * l1 * l2))) /
          2.;
  out2(3, 6) =
      -(copt1003 * copt479 * copt480 * copt481 * copt482 * copt732 *
        (2 * copt454 * copt459 * copt461 + 2 * copt304 * copt310 * copt475 -
         2 * copt471 * copt964)) /
          2. +
      (copt479 * copt480 * copt481 * copt482 * copt737 *
       (copt1439 * copt471 - 2 * copt459 * copt461 * copt730 +
        copt722 * copt964 -
        copt475 * (copt1432 * copt728 * l0 * l1 + copt1411 * copt726 * l0 * l2 +
                   copt1391 * copt724 * l1 * l2))) /
          2.;
  out2(3, 7) =
      -(copt1003 * copt479 * copt480 * copt481 * copt482 * copt732 *
        (2 * copt454 * copt459 * copt465 + 2 * copt304 * copt440 * copt475 -
         2 * copt471 * copt972)) /
          2. +
      (copt479 * copt480 * copt481 * copt482 * copt737 *
       (copt1507 * copt471 - 2 * copt459 * copt465 * copt730 +
        copt722 * copt972 -
        copt475 * (copt1500 * copt728 * l0 * l1 + copt1480 * copt726 * l0 * l2 +
                   copt1462 * copt724 * l1 * l2))) /
          2.;
  out2(3, 8) =
      -(copt1003 * copt479 * copt480 * copt481 * copt482 * copt732 *
        (2 * copt454 * copt459 * copt469 + 2 * copt304 * copt452 * copt475 -
         2 * copt471 * copt979)) /
          2. +
      (copt479 * copt480 * copt481 * copt482 * copt737 *
       (copt1576 * copt471 - 2 * copt459 * copt469 * copt730 +
        copt722 * copt979 -
        copt475 * (copt1569 * copt728 * l0 * l1 + copt1549 * copt726 * l0 * l2 +
                   copt1530 * copt724 * l1 * l2))) /
          2.;
  out2(3, 9) = (copt479 * copt480 * copt481 * copt482 * copt737 *
                (copt1586 * copt471 * copt581 * copt582 * l1 * l2 -
                 copt1586 * copt475 * copt724 * l1 * l2)) /
               2.;
  out2(3, 10) = (copt479 * copt480 * copt481 * copt482 * copt737 *
                 (copt1600 * copt471 * copt581 * copt582 * l1 * l2 -
                  copt1600 * copt475 * copt724 * l1 * l2)) /
                2.;
  out2(3, 11) = (copt479 * copt480 * copt481 * copt482 * copt737 *
                 (copt1611 * copt471 * copt581 * copt582 * l1 * l2 -
                  copt1611 * copt475 * copt724 * l1 * l2)) /
                2.;
  out2(3, 12) = (copt479 * copt480 * copt481 * copt482 * copt737 *
                 (copt1621 * copt471 * copt651 * copt652 * l0 * l2 -
                  copt1621 * copt475 * copt726 * l0 * l2)) /
                2.;
  out2(3, 13) = (copt479 * copt480 * copt481 * copt482 * copt737 *
                 (copt1631 * copt471 * copt651 * copt652 * l0 * l2 -
                  copt1631 * copt475 * copt726 * l0 * l2)) /
                2.;
  out2(3, 14) = (copt479 * copt480 * copt481 * copt482 * copt737 *
                 (copt1641 * copt471 * copt651 * copt652 * l0 * l2 -
                  copt1641 * copt475 * copt726 * l0 * l2)) /
                2.;
  out2(3, 15) = (copt479 * copt480 * copt481 * copt482 * copt737 *
                 (copt1651 * copt471 * copt719 * copt720 * l0 * l1 -
                  copt1651 * copt475 * copt728 * l0 * l1)) /
                2.;
  out2(3, 16) = (copt479 * copt480 * copt481 * copt482 * copt737 *
                 (copt1662 * copt471 * copt719 * copt720 * l0 * l1 -
                  copt1662 * copt475 * copt728 * l0 * l1)) /
                2.;
  out2(3, 17) = (copt479 * copt480 * copt481 * copt482 * copt737 *
                 (copt1673 * copt471 * copt719 * copt720 * l0 * l1 -
                  copt1673 * copt475 * copt728 * l0 * l1)) /
                2.;
  out2(4, 0) =
      (copt479 * copt480 * copt481 * copt482 * copt743 * copt807 *
       (-(copt814 * (2 * copt1702 * copt457 * copt459 + copt1713 * copt815 +
                     copt1700 * copt832)) -
        copt1069 * copt838 + copt1710 * copt849 + copt1723 * copt858)) /
          2. -
      (copt1694 * copt1696 * copt479 * copt480 * copt481 * copt482 * copt743 *
       copt860) /
          2.;
  out2(4, 1) =
      (copt479 * copt480 * copt481 * copt482 * copt743 * copt807 *
       (-(copt814 * (2 * copt1736 * copt457 * copt459 + copt1747 * copt815 +
                     copt1734 * copt832)) -
        copt1117 * copt838 + copt1744 * copt849 + copt1757 * copt858)) /
          2. -
      (copt1696 * copt1730 * copt479 * copt480 * copt481 * copt482 * copt743 *
       copt860) /
          2.;
  out2(4, 2) =
      (copt479 * copt480 * copt481 * copt482 * copt743 * copt807 *
       (-(copt814 * (2 * copt1771 * copt457 * copt459 + copt1782 * copt815 +
                     copt1769 * copt832)) -
        copt1165 * copt838 + copt1779 * copt849 + copt1792 * copt858)) /
          2. -
      (copt1696 * copt1765 * copt479 * copt480 * copt481 * copt482 * copt743 *
       copt860) /
          2.;
  out2(4, 3) =
      (copt479 * copt480 * copt481 * copt482 * copt743 * copt807 *
       (-(copt814 * (2 * copt308 * copt457 * copt459 + copt1821 * copt815)) -
        copt1234 * copt838 + copt1829 * copt849 + copt1835 * copt858)) /
          2. -
      (copt1696 * copt1818 * copt479 * copt480 * copt481 * copt482 * copt743 *
       copt860) /
          2.;
  out2(4, 4) =
      (copt479 * copt480 * copt481 * copt482 * copt743 * copt807 *
       (-(copt814 * (2 * copt438 * copt457 * copt459 + copt1865 * copt815)) -
        copt1299 * copt838 + copt1873 * copt849 + copt1879 * copt858)) /
          2. -
      (copt1696 * copt1862 * copt479 * copt480 * copt481 * copt482 * copt743 *
       copt860) /
          2.;
  out2(4, 5) =
      (copt479 * copt480 * copt481 * copt482 * copt743 * copt807 *
       (-(copt814 * (2 * copt450 * copt457 * copt459 + copt1905 * copt815)) -
        copt1369 * copt838 + copt1913 * copt849 + copt1919 * copt858)) /
          2. -
      (copt1696 * copt1902 * copt479 * copt480 * copt481 * copt482 * copt743 *
       copt860) /
          2.;
  out2(4, 6) =
      (copt479 * copt480 * copt481 * copt482 * copt743 * copt807 *
       (-(copt814 * (2 * copt301 * copt457 * copt459 + copt1947 * copt832)) -
        copt1439 * copt838 + copt1955 * copt849 + copt1961 * copt858)) /
          2. -
      (copt1696 * copt1944 * copt479 * copt480 * copt481 * copt482 * copt743 *
       copt860) /
          2.;
  out2(4, 7) =
      (copt479 * copt480 * copt481 * copt482 * copt743 * copt807 *
       (-(copt814 * (2 * copt416 * copt457 * copt459 + copt1988 * copt832)) -
        copt1507 * copt838 + copt1996 * copt849 + copt2002 * copt858)) /
          2. -
      (copt1696 * copt1985 * copt479 * copt480 * copt481 * copt482 * copt743 *
       copt860) /
          2.;
  out2(4, 8) =
      (copt479 * copt480 * copt481 * copt482 * copt743 * copt807 *
       (-(copt814 * (2 * copt447 * copt457 * copt459 + copt2027 * copt832)) -
        copt1576 * copt838 + copt2035 * copt849 + copt2041 * copt858)) /
          2. -
      (copt1696 * copt2024 * copt479 * copt480 * copt481 * copt482 * copt743 *
       copt860) /
          2.;
  out2(4, 9) = (copt479 * copt480 * copt481 * copt482 * copt743 * copt807 *
                (-(copt1586 * copt581 * copt582 * copt838 * l1 * l2) +
                 copt1586 * copt840 * copt858 * l1 * l2)) /
               2.;
  out2(4, 10) = (copt479 * copt480 * copt481 * copt482 * copt743 * copt807 *
                 (-(copt1600 * copt581 * copt582 * copt838 * l1 * l2) +
                  copt1600 * copt840 * copt858 * l1 * l2)) /
                2.;
  out2(4, 11) = (copt479 * copt480 * copt481 * copt482 * copt743 * copt807 *
                 (-(copt1611 * copt581 * copt582 * copt838 * l1 * l2) +
                  copt1611 * copt840 * copt858 * l1 * l2)) /
                2.;
  out2(4, 12) = (copt479 * copt480 * copt481 * copt482 * copt743 * copt807 *
                 (-(copt1621 * copt651 * copt652 * copt838 * l0 * l2) +
                  copt1621 * copt843 * copt858 * l0 * l2)) /
                2.;
  out2(4, 13) = (copt479 * copt480 * copt481 * copt482 * copt743 * copt807 *
                 (-(copt1631 * copt651 * copt652 * copt838 * l0 * l2) +
                  copt1631 * copt843 * copt858 * l0 * l2)) /
                2.;
  out2(4, 14) = (copt479 * copt480 * copt481 * copt482 * copt743 * copt807 *
                 (-(copt1641 * copt651 * copt652 * copt838 * l0 * l2) +
                  copt1641 * copt843 * copt858 * l0 * l2)) /
                2.;
  out2(4, 15) = (copt479 * copt480 * copt481 * copt482 * copt743 * copt807 *
                 (-(copt1651 * copt719 * copt720 * copt838 * l0 * l1) +
                  copt1651 * copt846 * copt858 * l0 * l1)) /
                2.;
  out2(4, 16) = (copt479 * copt480 * copt481 * copt482 * copt743 * copt807 *
                 (-(copt1662 * copt719 * copt720 * copt838 * l0 * l1) +
                  copt1662 * copt846 * copt858 * l0 * l1)) /
                2.;
  out2(4, 17) = (copt479 * copt480 * copt481 * copt482 * copt743 * copt807 *
                 (-(copt1673 * copt719 * copt720 * copt838 * l0 * l1) +
                  copt1673 * copt846 * copt858 * l0 * l1)) /
                2.;
  out2(5, 0) = (copt479 * copt480 * copt481 * copt482 * copt743 * copt807 *
                (copt1710 * copt814 + copt1069 * copt858 -
                 copt849 * (2 * copt1702 * copt304 * copt57 +
                            copt1713 * copt862 + copt1700 * copt865) -
                 copt1723 * copt867)) /
                   2. -
               (copt1694 * copt1696 * copt479 * copt480 * copt481 * copt482 *
                copt743 * copt870) /
                   2.;
  out2(5, 1) = (copt479 * copt480 * copt481 * copt482 * copt743 * copt807 *
                (copt1744 * copt814 + copt1117 * copt858 -
                 copt849 * (2 * copt1736 * copt304 * copt57 +
                            copt1747 * copt862 + copt1734 * copt865) -
                 copt1757 * copt867)) /
                   2. -
               (copt1696 * copt1730 * copt479 * copt480 * copt481 * copt482 *
                copt743 * copt870) /
                   2.;
  out2(5, 2) = (copt479 * copt480 * copt481 * copt482 * copt743 * copt807 *
                (copt1779 * copt814 + copt1165 * copt858 -
                 copt849 * (2 * copt1771 * copt304 * copt57 +
                            copt1782 * copt862 + copt1769 * copt865) -
                 copt1792 * copt867)) /
                   2. -
               (copt1696 * copt1765 * copt479 * copt480 * copt481 * copt482 *
                copt743 * copt870) /
                   2.;
  out2(5, 3) =
      (copt479 * copt480 * copt481 * copt482 * copt743 * copt807 *
       (copt1829 * copt814 + copt1234 * copt858 -
        copt849 * (2 * copt304 * copt308 * copt57 + copt1821 * copt862) -
        copt1835 * copt867)) /
          2. -
      (copt1696 * copt1818 * copt479 * copt480 * copt481 * copt482 * copt743 *
       copt870) /
          2.;
  out2(5, 4) =
      (copt479 * copt480 * copt481 * copt482 * copt743 * copt807 *
       (copt1873 * copt814 + copt1299 * copt858 -
        copt849 * (2 * copt304 * copt438 * copt57 + copt1865 * copt862) -
        copt1879 * copt867)) /
          2. -
      (copt1696 * copt1862 * copt479 * copt480 * copt481 * copt482 * copt743 *
       copt870) /
          2.;
  out2(5, 5) =
      (copt479 * copt480 * copt481 * copt482 * copt743 * copt807 *
       (copt1913 * copt814 + copt1369 * copt858 -
        copt849 * (2 * copt304 * copt450 * copt57 + copt1905 * copt862) -
        copt1919 * copt867)) /
          2. -
      (copt1696 * copt1902 * copt479 * copt480 * copt481 * copt482 * copt743 *
       copt870) /
          2.;
  out2(5, 6) =
      (copt479 * copt480 * copt481 * copt482 * copt743 * copt807 *
       (copt1955 * copt814 + copt1439 * copt858 -
        copt849 * (2 * copt301 * copt304 * copt57 + copt1947 * copt865) -
        copt1961 * copt867)) /
          2. -
      (copt1696 * copt1944 * copt479 * copt480 * copt481 * copt482 * copt743 *
       copt870) /
          2.;
  out2(5, 7) =
      (copt479 * copt480 * copt481 * copt482 * copt743 * copt807 *
       (copt1996 * copt814 + copt1507 * copt858 -
        copt849 * (2 * copt304 * copt416 * copt57 + copt1988 * copt865) -
        copt2002 * copt867)) /
          2. -
      (copt1696 * copt1985 * copt479 * copt480 * copt481 * copt482 * copt743 *
       copt870) /
          2.;
  out2(5, 8) =
      (copt479 * copt480 * copt481 * copt482 * copt743 * copt807 *
       (copt2035 * copt814 + copt1576 * copt858 -
        copt849 * (2 * copt304 * copt447 * copt57 + copt2027 * copt865) -
        copt2041 * copt867)) /
          2. -
      (copt1696 * copt2024 * copt479 * copt480 * copt481 * copt482 * copt743 *
       copt870) /
          2.;
  out2(5, 9) = (copt479 * copt480 * copt481 * copt482 * copt743 * copt807 *
                (copt1586 * copt581 * copt582 * copt858 * l1 * l2 -
                 copt1586 * copt840 * copt867 * l1 * l2)) /
               2.;
  out2(5, 10) = (copt479 * copt480 * copt481 * copt482 * copt743 * copt807 *
                 (copt1600 * copt581 * copt582 * copt858 * l1 * l2 -
                  copt1600 * copt840 * copt867 * l1 * l2)) /
                2.;
  out2(5, 11) = (copt479 * copt480 * copt481 * copt482 * copt743 * copt807 *
                 (copt1611 * copt581 * copt582 * copt858 * l1 * l2 -
                  copt1611 * copt840 * copt867 * l1 * l2)) /
                2.;
  out2(5, 12) = (copt479 * copt480 * copt481 * copt482 * copt743 * copt807 *
                 (copt1621 * copt651 * copt652 * copt858 * l0 * l2 -
                  copt1621 * copt843 * copt867 * l0 * l2)) /
                2.;
  out2(5, 13) = (copt479 * copt480 * copt481 * copt482 * copt743 * copt807 *
                 (copt1631 * copt651 * copt652 * copt858 * l0 * l2 -
                  copt1631 * copt843 * copt867 * l0 * l2)) /
                2.;
  out2(5, 14) = (copt479 * copt480 * copt481 * copt482 * copt743 * copt807 *
                 (copt1641 * copt651 * copt652 * copt858 * l0 * l2 -
                  copt1641 * copt843 * copt867 * l0 * l2)) /
                2.;
  out2(5, 15) = (copt479 * copt480 * copt481 * copt482 * copt743 * copt807 *
                 (copt1651 * copt719 * copt720 * copt858 * l0 * l1 -
                  copt1651 * copt846 * copt867 * l0 * l1)) /
                2.;
  out2(5, 16) = (copt479 * copt480 * copt481 * copt482 * copt743 * copt807 *
                 (copt1662 * copt719 * copt720 * copt858 * l0 * l1 -
                  copt1662 * copt846 * copt867 * l0 * l1)) /
                2.;
  out2(5, 17) = (copt479 * copt480 * copt481 * copt482 * copt743 * copt807 *
                 (copt1673 * copt719 * copt720 * copt858 * l0 * l1 -
                  copt1673 * copt846 * copt867 * l0 * l1)) /
                2.;
  return std::make_tuple(grad, val);
}

#endif  // hylc_strain_II
