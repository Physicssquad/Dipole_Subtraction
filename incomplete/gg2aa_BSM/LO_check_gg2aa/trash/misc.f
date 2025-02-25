      function dot(p,q)
      implicit double precision (a-h,o-z)
      dimension p(0:3),q(0:3)
      dot=p(0)*q(0)-p(1)*q(1)-p(2)*q(2)-p(3)*q(3)
      return
      end
