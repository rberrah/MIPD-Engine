(function () {
  let sender = null;
  let senderOrigin = '';
  let pending = '';
  let lastAck = null;
  function register() {
    if (!window.Shiny) return;
    window.Shiny.addCustomMessageHandler('workbench-ack', function (payload) {
      lastAck = { type: 'pk-workbench-ack', ...payload };
      if (sender) sender.postMessage(lastAck, senderOrigin);
      pending = '';
    });
    window.Shiny.addCustomMessageHandler('workbench-result', function (payload) {
      if (sender && lastAck?.id === payload.id) sender.postMessage({type: 'pk-workbench-result', ...payload}, senderOrigin);
    });
  }
  document.addEventListener('DOMContentLoaded', register);
  window.addEventListener('message', function (event) {
    const query = new URLSearchParams(window.location.search);
    if (query.get('bridge') !== 'workbench' || event.source !== window.opener || event.origin !== query.get('origin')) return;
    const origin = new URL(event.origin);
    const allowed = origin.origin === 'https://rberrah.github.io' ||
      (['localhost', '127.0.0.1'].includes(origin.hostname) && ['http:', 'https:'].includes(origin.protocol));
    const data = event.data;
    if (!allowed || !data || data.type !== 'pk-workbench' || typeof data.id !== 'string') return;
    if (!window.Shiny?.shinyapp?.config?.sessionId || window.Shiny.shinyapp.$socket?.readyState !== 1) return;
    if (lastAck?.id === data.id) { event.source.postMessage(lastAck, event.origin); return; }
    if (pending === data.id) return;
    sender = event.source;
    senderOrigin = event.origin;
    pending = data.id;
    register();
    window.Shiny.setInputValue('workbench_payload', data, { priority: 'event' });
  });
})();
