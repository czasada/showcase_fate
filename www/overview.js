// Card wheel of the Overview tab: arrows move by one card and wrap around at both ends.
function ovWheelStep(el) {
  var cards = el.children;
  if (cards.length < 2) return el.clientWidth;
  return cards[1].offsetLeft - cards[0].offsetLeft;
}

function ovWheel(id, dir) {
  var el = document.getElementById(id);
  if (!el) return;
  var max = el.scrollWidth - el.clientWidth;
  if (dir > 0 && el.scrollLeft >= max - 2) {
    el.scrollTo({ left: 0 });
  } else if (dir < 0 && el.scrollLeft <= 2) {
    el.scrollTo({ left: max });
  } else {
    el.scrollBy({ left: dir * ovWheelStep(el) });
  }
}

function ovWheelInit(id, posId) {
  var el = document.getElementById(id);
  var pos = document.getElementById(posId);
  if (!el || !pos) return;
  var update = function () {
    var step = ovWheelStep(el);
    var n = el.children.length;
    if (step <= 0 || n === 0) return;
    var first = Math.round(el.scrollLeft / step);
    var visible = Math.max(1, Math.round((el.clientWidth + 16) / step));
    pos.textContent = (first + 1) + "–" + Math.min(first + visible, n) + " of " + n;
  };
  el.addEventListener("scroll", function () { window.requestAnimationFrame(update); });
  window.addEventListener("resize", update);
  setTimeout(update, 50);
}
